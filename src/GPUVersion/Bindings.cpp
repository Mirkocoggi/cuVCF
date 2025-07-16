#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <cuda_fp16.h>
#include <pybind11/numpy.h>
#include <cstring>

#include "Parser.h"      // Contains vcf_parsed definition and run function (in Parser.cu)
#include "DataFrames.h"  // Contains classes: var_columns_df, alt_columns_df, sample_columns_df, alt_format_df
#include "DataStructures.h"

namespace py = pybind11;

/**
 * @class GTWrapper
 * @brief Wrapper class for genotype values with string-to-char mapping
 */
class GTWrapper {
    public:
        /** @brief Character value representing the genotype */
        char gt_value;
        
        /** @brief Shared map for string-to-char conversion */
        static std::map<std::string, char> GTMap;
        
        GTWrapper(char value) : gt_value(value) {}
        
        /**
         * @brief Converts the stored genotype value back to its string representation
         * @return String representation of the genotype
         */
        std::string to_string() const {
            for (const auto& pair : GTMap) {
                if (pair.second == gt_value) {
                    return pair.first;
                }
            }
            return std::string(1, gt_value);
        }
    };
    
    std::map<std::string, char> GTWrapper::GTMap;
    
    void init_GTMap() {
        int value = 0;
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key1 = std::to_string(i) + "|" + std::to_string(j);
                GTWrapper::GTMap[key1] = value;
                value++;
            }
        }
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key2 = std::to_string(i) + "/" + std::to_string(j);
                GTWrapper::GTMap[key2] = value;
                value++;
            }
        }
        GTWrapper::GTMap[".|."] = static_cast<char>(254);
        GTWrapper::GTMap["./."] = static_cast<char>(255);
    }
    
    void bind_GTWrapper(py::module &m) {
        py::class_<GTWrapper>(m, "GT")
            .def(py::init<char>())
            .def("__repr__", &GTWrapper::to_string);
    }
    

    struct half_wrapper {
        __half value;
        half_wrapper() : value(__float2half(0.0f)) {}
        half_wrapper(float f) : value(__float2half(f)) {}
        half_wrapper(const __half &h) : value(h) {}
        float to_float() const { return __half2float(value); }
    };
    
    PYBIND11_MAKE_OPAQUE(std::vector<half_wrapper>);
    
    void bind_vector_half(py::module &m) {
        py::bind_vector<std::vector<half_wrapper>>(m, "VectorHalf");
    }
    
    std::vector<half_wrapper> convert_half_vector(const std::vector<__half>& src) {
        std::vector<half_wrapper> result;
        result.reserve(src.size());
        for (const auto &h : src) {
            result.push_back(half_wrapper(h));
        }
        return result;
    }
    
    void set_half_vector(std::vector<__half>& dest, const std::vector<half_wrapper>& src) {
        dest.clear();
        for (const auto &hw : src) {
            dest.push_back(hw.value);
        }
    }
    
/**
 * @brief Converts a vector of CUDA half values to a NumPy array of uint16
 * 
 * @param vec Vector of CUDA half-precision values
 * @return py::array_t<uint16_t> NumPy array containing the raw bits of half values
 * 
 * @details Creates a memory-managed NumPy array from half-precision values,
 * preserving the bit representation for Python-side processing
 */
py::array_t<uint16_t> half_vector_to_uint16(const std::vector<__half>& vec) {
    // Allocate heap storage for raw data
    auto* raw = new std::vector<uint16_t>();
    raw->reserve(vec.size());

    // Convert each half value to its raw bits
    for (const __half &h : vec) {
        uint16_t bits = *reinterpret_cast<const uint16_t*>(&h);
        raw->push_back(bits);
    }

    // Create capsule for memory management
    py::capsule free_when_done(raw, [](void *f) {
        auto* ptr = reinterpret_cast<std::vector<uint16_t>*>(f);
        delete ptr;
    });

    // Create and return NumPy array
    return py::array_t<uint16_t>(
        { raw->size() },            // shape
        { sizeof(uint16_t) },       // strides
        raw->data(),                // data pointer
        free_when_done              // memory management
    );
}

py::dict get_var_columns_data(const var_columns_df &df) {
    py::dict d;

    // var_number (unsigned int vector)
    d["var_number"] = py::array_t<unsigned int>(df.var_number.size(), df.var_number.data());
    
    d["chrom"] = py::array_t<char>(df.chrom.size(), df.chrom.data());
    
    d["pos"] = py::array_t<unsigned int>(df.pos.size(), df.pos.data());
    
    d["id"] = py::cast(df.id);
    d["ref"] = py::cast(df.ref);
    d["filter"] = py::cast(df.filter);
    
    d["qual"] = half_vector_to_uint16(df.qual);

    for (const auto &info : df.in_float) {
        d[info.name.c_str()] = half_vector_to_uint16(info.i_float);
    }

    for (const auto &info : df.in_flag) {
        if (std::any_of(info.i_flag.begin(), info.i_flag.end(), [](auto val) { return val == 1; })) {
            d[info.name.c_str()] = py::cast(info.i_flag);
        }
    }
    for (const auto &info : df.in_string) {
        d[info.name.c_str()] = py::cast(info.i_string);
    }
    for (const auto &info : df.in_int) {
        d[info.name.c_str()] = py::cast(info.i_int);
    }

    return d;
}


py::dict get_alt_columns_data(const alt_columns_df &df) {
    py::dict d;

    d["var_id"] = py::array_t<unsigned int>(df.var_id.size(), df.var_id.data());
    d["alt_id"] = py::array_t<unsigned char>(df.alt_id.size(), df.alt_id.data());
    d["alt"] = py::cast(df.alt);

    for (const auto &info : df.alt_float) {
        d[info.name.c_str()] = half_vector_to_uint16(info.i_float);
    }
    
    for (const auto &info : df.alt_flag) {
        if (std::any_of(info.i_flag.begin(), info.i_flag.end(), [](auto val) { return val == 1; })) {
            d[info.name.c_str()] = py::cast(info.i_flag);
        }
    }

    for (const auto &info : df.alt_string) {
        d[info.name.c_str()] = py::cast(info.i_string);
    }
    for (const auto &info : df.alt_int) {
        d[info.name.c_str()] = py::cast(info.i_int);
    }

    return d;
}

py::dict get_sample_columns_data(const sample_columns_df &df) {
    py::dict d;
    
    d["var_id"] = py::array_t<unsigned int>(df.var_id.size(), df.var_id.data());

    d["samp_id"] = py::array_t<unsigned short>(df.samp_id.size(), df.samp_id.data());
    
    for (const auto &s : df.samp_float) {
        if (!s.i_float.empty()) {
            d[s.name.c_str()] = half_vector_to_uint16(s.i_float);
        }
    }
    
    for (const auto &s : df.samp_flag) {
        if (!s.i_flag.empty() &&
            std::any_of(s.i_flag.begin(), s.i_flag.end(), [](auto val) { return val == 1; })) {
            d[s.name.c_str()] = py::cast(s.i_flag);
        }
    }

    for (const auto &s : df.samp_string) {
        if (!s.i_string.empty()) {
            d[s.name.c_str()] = py::cast(s.i_string);
        }
    }
    for (const auto &s : df.samp_int) {
        if (!s.i_int.empty()) {
            d[s.name.c_str()] = py::cast(s.i_int);
        }
    }
    
    int i = 0;
    for (const auto &gt_struct : df.sample_GT) {
        if (!gt_struct.GT.empty()) {
            std::string key = "GT" + std::to_string(i);
            d[key.c_str()] = py::cast(gt_struct.GT);
            i++;
        }
    }
    
    return d;
}

py::dict get_alt_format_data(const alt_format_df &df) {
    py::dict d;
    
    d["var_id"] = py::array_t<unsigned int>(df.var_id.size(), df.var_id.data());
    d["samp_id"] = py::array_t<unsigned short>(df.samp_id.size(), df.samp_id.data());
    d["alt_id"] = py::array_t<char>(df.alt_id.size(), df.alt_id.data());
    
    for (const auto &s : df.samp_float) {
        if (!s.i_float.empty()) {
            d[s.name.c_str()] = half_vector_to_uint16(s.i_float);
        }
    }

    for (const auto &s : df.samp_flag) {
        if (!s.i_flag.empty() &&
            std::any_of(s.i_flag.begin(), s.i_flag.end(), [](auto val) { return val == 1; })) {
            d[s.name.c_str()] = py::cast(s.i_flag);
        }
    }
    for (const auto &s : df.samp_string) {
        if (!s.i_string.empty()) {
            d[s.name.c_str()] = py::cast(s.i_string);
        }
    }
    for (const auto &s : df.samp_int) {
        if (!s.i_int.empty()) {
            d[s.name.c_str()] = py::cast(s.i_int);
        }
    }
    
    if (!df.sample_GT.GT.empty()) {
        d["GT"] = py::cast(df.sample_GT.GT);
    }
    
    return d;
}

/**
 * @brief Binds VCF data structures to Python
 * 
 * Creates Python bindings for the following C++ classes:
 * - var_columns_df: Core variant data
 * - alt_columns_df: Alternative allele data
 * - sample_columns_df: Sample-specific data
 * - alt_format_df: Alternative format data
 * 
 * Each class is exposed with its members and methods to Python,
 * with proper type conversion between C++ and Python types.
 */
PYBIND11_MODULE(GPUParser, m) {
    m.doc() = "Python bindings for CUDA-accelerated VCF parser using pybind11";
    py::bind_map<std::map<std::string, char>>(m, "FilterMap");

    py::class_<half_wrapper>(m, "half")
        .def(py::init<>())
        .def(py::init<float>())
        .def("__float__", &half_wrapper::to_float)
        .def("__repr__", [](const half_wrapper &h) {
            return std::to_string(h.to_float());
        });
    bind_vector_half(m);

    init_GTMap();
    bind_GTWrapper(m);
    m.attr("GTMapGlobal") = py::cast(&GTWrapper::GTMap);

    py::class_<samp_GT>(m, "samp_GT")
        .def(py::init<>())
        .def_readwrite("GT", &samp_GT::GT)
        .def_readwrite("numb", &samp_GT::numb); 

    py::class_<info_flag>(m, "info_flag")
        .def(py::init<>())
        .def_readwrite("i_flag", &info_flag::i_flag)
        .def_readwrite("name", &info_flag::name);

    py::class_<info_string>(m, "info_string")
        .def(py::init<>())
        .def_readwrite("i_string", &info_string::i_string)
        .def_readwrite("name", &info_string::name);

    py::class_<info_float>(m, "info_float")
        .def(py::init<>())
        .def_readwrite("i_float", &info_float::i_float)
        .def_readwrite("name", &info_float::name);

    py::class_<info_int>(m, "info_int")
        .def(py::init<>())
        .def_readwrite("i_int", &info_int::i_int)
        .def_readwrite("name", &info_int::name);

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

    py::class_<var_columns_df>(m, "var_columns_df")
        .def(py::init<>())
        .def_readwrite("var_number", &var_columns_df::var_number)
        .def_readwrite("chrom", &var_columns_df::chrom)
        .def_readwrite("pos", &var_columns_df::pos)
        .def_readwrite("id", &var_columns_df::id)
        .def_readwrite("ref", &var_columns_df::ref)
        .def_property("qual",
            [](const var_columns_df &self) -> std::vector<half_wrapper> {
                return convert_half_vector(self.qual);
            },
            [](var_columns_df &self, const std::vector<half_wrapper>& vec) {
                set_half_vector(self.qual, vec);
            },
            "Quality scores as half values wrapped in half_wrapper")
        .def_readwrite("filter", &var_columns_df::filter)
        .def_readwrite("filter_map", &var_columns_df::filter_map)
        .def_readwrite("in_float", &var_columns_df::in_float)
        .def_readwrite("in_flag", &var_columns_df::in_flag)
        .def_readwrite("in_string", &var_columns_df::in_string)
        .def_readwrite("in_int", &var_columns_df::in_int)
        .def_readwrite("info_map1", &var_columns_df::info_map1)
        .def("print", &var_columns_df::print, py::arg("num_lines"), "Print variant columns");

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
        .def("print", &alt_columns_df::print, py::arg("n"), "Print alternative columns");

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
        .def("print", &sample_columns_df::print, py::arg("n"), "Print sample columns");

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
        .def("print", &alt_format_df::print, py::arg("n"), "Print alternative formatted columns");

    py::class_<vcf_parsed>(m, "vcf_parsed")
        .def(py::init<>())
        .def("run", [](vcf_parsed &self, const std::string &vcf_filename, int num_threadss) {
            char* c_vcf_filename = new char[vcf_filename.size() + 1];
            std::strcpy(c_vcf_filename, vcf_filename.c_str());
            self.run(c_vcf_filename, num_threadss);
            delete[] c_vcf_filename;
        }, py::arg("vcf_filename"), py::arg("num_threadss"),
           "Parses the VCF file using the specified number of threads")
        .def_readwrite("id", &vcf_parsed::id)
        .def_readwrite("header", &vcf_parsed::header)
        .def_readwrite("INFO", &vcf_parsed::INFO)
        .def_readwrite("FORMAT", &vcf_parsed::FORMAT)
        .def_readwrite("var_columns", &vcf_parsed::var_columns)
        .def_readwrite("alt_columns", &vcf_parsed::alt_columns)
        .def_readwrite("samp_columns", &vcf_parsed::samp_columns)
        .def_readwrite("alt_sample", &vcf_parsed::alt_sample)
        .def("print_header", &vcf_parsed::print_header);

        m.def("get_var_columns_data", &get_var_columns_data, 
            "Returns a dict of NumPy arrays containing var_columns_df structure data");

        m.def("get_alt_columns_data", &get_alt_columns_data, 
            "Returns a dict of NumPy arrays containing alt_columns_df structure data");

        m.def("get_sample_columns_data", &get_sample_columns_data,
            "Returns a dict of NumPy arrays containing sample_columns_df structure data");

        m.def("get_alt_format_data", &get_alt_format_data,
            "Returns a dict of NumPy arrays containing alt_format_df structure data");

}
