/**
 * @file CPUParser.cpp
 * @brief CPU‑only pybind11 bindings for the multi‑threaded VCF parser
 *        and its supporting data structures.
 *
 * This translation unit exposes to Python a complete, header‑only, CPU
 * implementation of the VCF parsing pipeline that is API‑compatible with
 * the CPU version (``GPUParser``).  In particular it provides:
 *   - A compact ``GTWrapper`` for genotype encoding/decoding.
 *   - A ``half_wrapper`` helper so half‑precision values can be exchanged
 *     with NumPy using plain ``uint16_t`` buffers.
 *   - Utility functions that serialise the column‑oriented C++ structures
 *     (``var_columns_df``, ``alt_columns_df`` …) into NumPy‑backed
 *     ``dict`` objects.
 *   - Full pybind11 bindings for every public structure required by the
 *     Python front‑end, including the top‑level ``vcf_parsed`` parser.
 *
 * @note All code paths are CPU‑only; no CUDA headers or kernels are
 *       referenced here.  The compilation unit can therefore be used on
 *       systems without an NVIDIA GPU or a CUDA toolchain installed.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <half.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>        
#include <pybind11/stl_bind.h>   
#include "VCF_parsed.h"
#include "VCF_var_columns_df.h"
#include "VCFparser_mt_col_struct.h"   

namespace py = pybind11;

/**
 * @class GTWrapper
 * @brief Lightweight value‑type that stores a genotype code.
 *
 * Genotypes such as "0|1" or "1/1" are encoded as a single ``char`` so they
 * can be stored compactly inside column‑wise buffers.  The static map
 * #GTMap is initialised by ::init_GTMap() at module initialisation time and
 * provides bidirectional mapping between the human‑readable string and the
 * compact numerical code.
 */
class GTWrapper {
  public:
    char gt_value;

    /**
     * Global lookup table mapping genotype strings to their numeric code.
     *
     * @warning Filled **once** by ::init_GTMap().  Do **not** modify after
     *          initialisation as the numeric values are assumed to remain
     *          stable for the lifetime of the program.
     */
    static std::map<std::string,char> GTMap;

    /**
     * @brief Construct from a raw numeric code.
     * @param v Encoded genotype value.
     */
    explicit GTWrapper(char v) : gt_value(v) {}

    /**
     * @brief Convert the stored code back to its string representation.
     * @return The canonical genotype string (e.g. "0|1").
     */
    std::string to_string() const {
        for (auto const& p : GTMap)
            if (p.second == gt_value) return p.first;
        return std::string(1, gt_value);
    }
};

std::map<std::string,char> GTWrapper::GTMap;

/**
 * @brief Populate ::GTWrapper::GTMap with every possible genotype up to 10/10.
 *
 * The mapping is bidirectional: for each genotype string the corresponding
 * compact ``char`` code is stored.  Codes ``254`` and ``255`` are reserved
 * for missing values (".|." and "./.").
 */
static void init_GTMap() {
    int v = 0;
    for (int i=0;i<11;++i)
        for (int j=0;j<11;++j)
            GTWrapper::GTMap[std::to_string(i)+'|'+std::to_string(j)] = v++;
    for (int i=0;i<11;++i)
        for (int j=0;j<11;++j)
            GTWrapper::GTMap[std::to_string(i)+'/'+std::to_string(j)] = v++;
    GTWrapper::GTMap[".|."] = static_cast<char>(254);
    GTWrapper::GTMap["./."] = static_cast<char>(255);
}

/**
 * @struct half_wrapper
 * @brief Tiny helper around ``half`` to satisfy pybind11's type requirements.
 *
 * pybind11 cannot bind ``half`` directly, but it can bind a POD struct.
 */
struct half_wrapper {
    half value;
    half_wrapper() : value(half(0.0f)) {}
    explicit half_wrapper(float f) : value(half(f)) {}
    explicit half_wrapper(const half& h): value(h) {}
    float to_float() const { return static_cast<float>(value); }
};

PYBIND11_MAKE_OPAQUE(std::vector<half_wrapper>);

/**
 * @brief Convert a ``std::vector<half>`` into a vector of @ref half_wrapper.
 *
 * @param src Source vector of IEEE‑754 half‑precision values.
 * @return Destination vector of wrapper objects, each wrapping the
 *         corresponding element of @p src.
 */
static std::vector<half_wrapper> convert_half_vector(const std::vector<half>& src){
    std::vector<half_wrapper> out; out.reserve(src.size());
    for (auto const& h:src) out.emplace_back(h);
    return out;
}

/**
 * @brief Copy data from a vector of @ref half_wrapper back to a vector of
 *        plain ``half``.
 *
 * @param[out] dst Destination vector (cleared then filled).
 * @param[in]  src Source vector containing wrapper objects.
 */
static void set_half_vector(std::vector<half>& dst,
                            const std::vector<half_wrapper>& src){
    dst.clear(); dst.reserve(src.size());
    for (auto const& w:src) dst.push_back(w.value);
}

/**
 * @brief Re‑interpret a ``std::vector<half>`` as a NumPy ``uint16`` array.
 *
 * The raw bits of each ``half`` are copied into an owning ``std::vector`` of
 * ``uint16_t`` which is then exposed to Python via a NumPy view.  A
 * ``py::capsule`` is attached to manage the lifetime of the vector.
 *
 * @param v Input vector of half‑precision numbers.
 * @return NumPy array of shape ``(v.size(),)`` and dtype ``uint16``.
 */
static py::array_t<uint16_t> half_vector_to_uint16(const std::vector<half>& v){
    const py::ssize_t n = static_cast<py::ssize_t>(v.size());

    if (n == 0) {
        return py::array_t<uint16_t>(0);
    }

    auto* raw = new std::vector<uint16_t>;
    raw->reserve(n);
    for (auto const& h : v)
        raw->push_back(*reinterpret_cast<const uint16_t*>(&h));

    py::capsule cap(raw, [](void* p){
        delete reinterpret_cast<std::vector<uint16_t>*>(p);
    });

    return py::array_t<uint16_t>(n, raw->data(), cap);                              // capsule = owner
}

/**
 * @brief Serialise a @ref var_columns_df into a Python ``dict`` of NumPy arrays.
 *
 * The function iterates over every sub‑field in the dataframe and converts it
 * into a suitable NumPy view or Python object.  Variadic INFO vectors are
 * added only if at least one element is non‑zero (for flags) or non‑empty.
 */
static py::dict get_var_columns_data(const var_columns_df& df){
    py::dict d;
    d["var_number"] = py::array_t<unsigned int>(df.var_number.size(),df.var_number.data());
    d["chrom"]      = py::array_t<char>        (df.chrom.size(),df.chrom.data());
    d["pos"]        = py::array_t<unsigned int>(df.pos.size(),df.pos.data());
    d["id"]         = py::cast(df.id);
    d["ref"]        = py::cast(df.ref);
    d["filter"]     = py::cast(df.filter);
    d["qual"]       = half_vector_to_uint16(df.qual);
    for (auto const& x: df.in_float)   d[x.name.c_str()] = half_vector_to_uint16(x.i_float);
    for (auto const& x: df.in_flag)
        if (std::any_of(x.i_flag.begin(),x.i_flag.end(),[](auto v){return v==1;}))
            d[x.name.c_str()] = py::cast(x.i_flag);
    for (auto const& x: df.in_string)  d[x.name.c_str()] = py::cast(x.i_string);
    for (auto const& x: df.in_int)     d[x.name.c_str()] = py::cast(x.i_int);
    return d;
}

/**
 * @brief Serialise an @ref alt_columns_df into a Python ``dict`` of NumPy arrays.
 */
static py::dict get_alt_columns_data(const alt_columns_df& df){
    py::dict d;
    d["var_id"] = py::array_t<unsigned int>(df.var_id.size(),df.var_id.data());
    d["alt_id"] = py::array_t<char>(df.alt_id.size(), df.alt_id.data());
    d["alt"]    = py::cast(df.alt);
    for (auto const& x: df.alt_float)  d[x.name.c_str()] = half_vector_to_uint16(x.i_float);
    for (auto const& x: df.alt_flag)
        if (std::any_of(x.i_flag.begin(),x.i_flag.end(),[](auto v){return v==1;}))
            d[x.name.c_str()] = py::cast(x.i_flag);
    for (auto const& x: df.alt_string) d[x.name.c_str()] = py::cast(x.i_string);
    for (auto const& x: df.alt_int)    d[x.name.c_str()] = py::cast(x.i_int);
    return d;
}

/**
 * @brief Serialise a @ref sample_columns_df into a Python ``dict`` of NumPy arrays.
 */
static py::dict get_sample_columns_data(const sample_columns_df& df){
    py::dict d;
    d["var_id"] = py::array_t<unsigned int>(df.var_id.size(),df.var_id.data());
    d["samp_id"]= py::array_t<unsigned short>(df.samp_id.size(),df.samp_id.data());
    for (auto const& s: df.samp_float)
        if (!s.i_float.empty()) d[s.name.c_str()] = half_vector_to_uint16(s.i_float);
    for (auto const& s: df.samp_flag)
        if (std::any_of(s.i_flag.begin(),s.i_flag.end(),[](auto v){return v==1;}))
            d[s.name.c_str()] = py::cast(s.i_flag);
    for (auto const& s: df.samp_string) if(!s.i_string.empty()) d[s.name.c_str()] = py::cast(s.i_string);
    for (auto const& s: df.samp_int)    if(!s.i_int.empty())    d[s.name.c_str()] = py::cast(s.i_int);
    int i=0;
    for (auto const& gt: df.sample_GT)
        if(!gt.GT.empty()) d[("GT"+std::to_string(i++)).c_str()] = py::cast(gt.GT);
    return d;
}

/**
 * @brief Serialise an @ref alt_format_df into a Python ``dict`` of NumPy arrays.
 */
static py::dict get_alt_format_data(const alt_format_df& df){
    py::dict d;
    d["var_id"] = py::array_t<unsigned int>(df.var_id.size(),df.var_id.data());
    d["samp_id"]= py::array_t<unsigned short>(df.samp_id.size(),df.samp_id.data());
    d["alt_id"] = py::array_t<char>(df.alt_id.size(),df.alt_id.data());
    for (auto const& s: df.samp_float)
        if (!s.i_float.empty()) d[s.name.c_str()] = half_vector_to_uint16(s.i_float);
    for (auto const& s: df.samp_flag)
        if (std::any_of(s.i_flag.begin(),s.i_flag.end(),[](auto v){return v==1;}))
            d[s.name.c_str()] = py::cast(s.i_flag);
    for (auto const& s: df.samp_string) if(!s.i_string.empty()) d[s.name.c_str()] = py::cast(s.i_string);
    for (auto const& s: df.samp_int)    if(!s.i_int.empty())    d[s.name.c_str()] = py::cast(s.i_int);
    if(!df.sample_GT.GT.empty()) d["GT"] = py::cast(df.sample_GT.GT);
    return d;
}

/**
 * @brief Top‑level pybind11 module initialisation function.
 *
 * Name: ``CPUParser`` (mirrors ``GPUParser``).  All classes and helper
 * functions defined above are bound here.  The module is intentionally kept
 * API‑compatible with its GPU counterpart so that the same Python code can
 * switch between CPU and GPU back‑ends by importing the appropriate module.
 */
PYBIND11_MODULE(CPUParser, m) {
    m.doc() = "CPU-only bindings – API compatibile con GPUParser";

    /*  basic types  */
    py::class_<half_wrapper>(m,"half")
        .def(py::init<>())
        .def(py::init<float>())
        .def("__float__",&half_wrapper::to_float)
        .def("__repr__",[](const half_wrapper& h){return std::to_string(h.to_float());});
    py::bind_vector<std::vector<half_wrapper>>(m,"VectorHalf");

    /*  GT wrapper  */
    init_GTMap();
    py::class_<GTWrapper>(m,"GT")
        .def(py::init<char>())
        .def("__repr__",&GTWrapper::to_string);
    py::bind_map<std::map<std::string,char>>(m,"FilterMap");
    m.attr("GTMapGlobal") = py::cast(&GTWrapper::GTMap);

    /*  helper NumPy dict  */
    m.def("get_var_columns_data",   &get_var_columns_data);
    m.def("get_alt_columns_data",   &get_alt_columns_data);
    m.def("get_sample_columns_data",&get_sample_columns_data);
    m.def("get_alt_format_data",    &get_alt_format_data);

    /*  --- bind all the structs (info_flag, samp_*, header_element …) ---  */
    py::class_<info_flag>(m,"info_flag").def(py::init<>())
        .def_readwrite("i_flag",&info_flag::i_flag).def_readwrite("name",&info_flag::name);
    py::class_<info_string>(m,"info_string").def(py::init<>())
        .def_readwrite("i_string",&info_string::i_string).def_readwrite("name",&info_string::name);
    py::class_<info_float>(m,"info_float").def(py::init<>())
        .def_readwrite("i_float",&info_float::i_float).def_readwrite("name",&info_float::name);
    py::class_<info_int>(m,"info_int").def(py::init<>())
        .def_readwrite("i_int",&info_int::i_int).def_readwrite("name",&info_int::name);

    py::class_<samp_Flag>(m, "samp_Flag")
        .def(py::init<>())
        .def_readwrite("i_flag", &samp_Flag::i_flag)
        .def_readwrite("name", &samp_Flag::name)
        .def_readwrite("numb", &samp_Flag::numb);

    py::class_<samp_String>(m, "samp_String")
        .def(py::init<>())
        .def_readwrite("i_string", &samp_String::i_string)
        .def_readwrite("name", &samp_String::name)
        .def_readwrite("numb", &samp_String::numb);

    py::class_<samp_Float>(m, "samp_Float")
        .def(py::init<>())
        .def_readwrite("i_float", &samp_Float::i_float)
        .def_readwrite("name", &samp_Float::name)
        .def_readwrite("numb", &samp_Float::numb);

    py::class_<samp_Int>(m, "samp_Int")
        .def(py::init<>())
        .def_readwrite("i_int", &samp_Int::i_int)
        .def_readwrite("name", &samp_Int::name)
        .def_readwrite("numb", &samp_Int::numb);

    py::class_<header_element>(m, "header_element")
        .def(py::init<>())
        .def_readwrite("ID", &header_element::ID)
        .def_readwrite("Number", &header_element::Number)
        .def_readwrite("Type", &header_element::Type)
        .def_readwrite("total_values", &header_element::total_values)
        .def_readwrite("alt_values", &header_element::alt_values)
        .def_readwrite("no_alt_values", &header_element::no_alt_values)
        .def_readwrite("ints_alt", &header_element::ints_alt)
        .def_readwrite("floats_alt", &header_element::floats_alt)
        .def_readwrite("strings_alt", &header_element::strings_alt)
        .def_readwrite("flags_alt", &header_element::flags_alt)
        .def_readwrite("ints", &header_element::ints)
        .def_readwrite("floats", &header_element::floats)
        .def_readwrite("strings", &header_element::strings)
        .def_readwrite("flags", &header_element::flags);

    py::class_<alt_columns_df>(m, "alt_columns_df")
        .def(py::init<>())
        .def_readwrite("var_id", &alt_columns_df::var_id)
        .def_readwrite("alt_id", &alt_columns_df::alt_id)
        .def_readwrite("alt", &alt_columns_df::alt)
        .def_readwrite("alt_float", &alt_columns_df::alt_float)
        .def_readwrite("alt_flag", &alt_columns_df::alt_flag)
        .def_readwrite("alt_string", &alt_columns_df::alt_string)
        .def_readwrite("alt_int", &alt_columns_df::alt_int)
        .def_readwrite("numAlt", &alt_columns_df::numAlt)
        .def("print", &alt_columns_df::print);

     py::class_<sample_columns_df>(m, "sample_columns_df")
        .def(py::init<>())
        .def_readwrite("var_id", &sample_columns_df::var_id)
        .def_readwrite("samp_id", &sample_columns_df::samp_id)
        .def_readwrite("samp_float", &sample_columns_df::samp_float)
        .def_readwrite("samp_flag", &sample_columns_df::samp_flag)
        .def_readwrite("samp_string", &sample_columns_df::samp_string)
        .def_readwrite("samp_int", &sample_columns_df::samp_int)
        .def_readwrite("sampNames", &sample_columns_df::sampNames)
        .def_readwrite("numSample", &sample_columns_df::numSample)
        .def_readwrite("sample_GT", &sample_columns_df::sample_GT)
        .def("print", &sample_columns_df::print);

    py::class_<alt_format_df>(m, "alt_format_df")
        .def(py::init<>())
        .def_readwrite("var_id", &alt_format_df::var_id)
        .def_readwrite("samp_id", &alt_format_df::samp_id)
        .def_readwrite("alt_id", &alt_format_df::alt_id)
        .def_readwrite("samp_float", &alt_format_df::samp_float)
        .def_readwrite("samp_flag", &alt_format_df::samp_flag)
        .def_readwrite("samp_string", &alt_format_df::samp_string)
        .def_readwrite("samp_int", &alt_format_df::samp_int)
        .def_readwrite("sampNames", &alt_format_df::sampNames)
        .def_readwrite("numSample", &alt_format_df::numSample)
        .def_readwrite("sample_GT", &alt_format_df::sample_GT)
        .def("print", &alt_format_df::print);

    /*  var_columns with getter/setter for qual  */
    py::class_<var_columns_df>(m,"var_columns_df")
        .def(py::init<>())
        .def_readwrite("var_number",&var_columns_df::var_number)
        .def_readwrite("chrom",&var_columns_df::chrom)
        .def_readwrite("pos",&var_columns_df::pos)
        .def_readwrite("id",&var_columns_df::id)
        .def_readwrite("ref",&var_columns_df::ref)
        .def_property("qual",
            [](const var_columns_df& self){ return convert_half_vector(self.qual); },
            [](var_columns_df& self,const std::vector<half_wrapper>& v){ set_half_vector(self.qual,v); },
            "quality (half-wrapper)")
        .def_readwrite("filter",&var_columns_df::filter)
        .def_readwrite("in_float",&var_columns_df::in_float)
        .def_readwrite("in_flag",&var_columns_df::in_flag)
        .def_readwrite("in_string",&var_columns_df::in_string)
        .def_readwrite("in_int",&var_columns_df::in_int)
        .def_readwrite("info_map1",&var_columns_df::info_map1)
        .def("print",&var_columns_df::print);

    /*  parser main object  */
    py::class_<vcf_parsed>(m,"vcf_parsed")
        .def(py::init<>())
        .def("run",&vcf_parsed::run,py::arg("vcf_filename"),py::arg("num_threads"))
        .def_readwrite("id",&vcf_parsed::id)
        .def_readwrite("header",&vcf_parsed::header)
        .def_readwrite("INFO",&vcf_parsed::INFO)
        .def_readwrite("FORMAT",&vcf_parsed::FORMAT)
        .def_readwrite("var_columns",&vcf_parsed::var_columns)
        .def_readwrite("alt_columns",&vcf_parsed::alt_columns)
        .def_readwrite("samp_columns",&vcf_parsed::samp_columns)
        .def_readwrite("alt_sample",&vcf_parsed::alt_sample)
        .def("print_header",&vcf_parsed::print_header);
}