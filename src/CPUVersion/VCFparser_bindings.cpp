#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <half.h>
#include "VCF_parsed.h"
#include "VCF_var_columns_df.h"
#include "VCFparser_mt_col_struct.h"
namespace py = pybind11;

class GTWrapper {
public:
    char gt_value;  // Il valore `char` rappresentativo in C++
    static std::map<std::string, char> GTMap;  // La mappa condivisa per la conversione

    // Costruttore
    GTWrapper(char value) : gt_value(value) {}

    // Funzione per rappresentare il valore `gt_value` come stringa utilizzando `GTMap`
    std::string to_string() const {
        for (const auto& pair : GTMap) {
            if (pair.second == gt_value) {
                return pair.first;  // Restituisce la stringa corrispondente a `char`
            }
        }
        return std::string(1, gt_value);  // Se non trovato, ritorna semplicemente `char`
    }
};

// Inizializza GTMap
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
    GTMap[".|."] = static_cast<char>(254);
    GTMap["./."] = static_cast<char>(255);
}

void bind_GTWrapper(py::module &m) {
    py::class_<GTWrapper>(m, "GT")
        .def(py::init<char>())  // Costruttore che accetta un `char`
        .def("__repr__", &GTWrapper::to_string);  // Rappresentazione leggibile usando `GTMap`
}


// Bind `half` as a Python type
void bind_half(py::module &m) {
    py::class_<half>(m, "half")
        .def(py::init<>())
        .def("__float__", [](const half &h) { return static_cast<float>(h); })  // Conversion to float
        .def("__repr__", [](const half &h) {
            return std::to_string(static_cast<float>(h));
        });
}

// Explicitly bind `std::vector<half>`
PYBIND11_MAKE_OPAQUE(std::vector<half>);

void bind_vector_half(py::module &m) {
    py::bind_vector<std::vector<half>>(m, "VectorHalf");
}

PYBIND11_MODULE(VCFparser_mt_col_struct, m) {
    bind_half(m);            // Register `half`
    bind_vector_half(m);     // Register `std::vecto
    init_GTMap();            // Inizializza GTMap
    bind_GTWrapper(m);       // Registra GTWrapper

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

    py::class_<var_columns_df>(m, "var_columns_df")
        .def(py::init<>())
        .def_readwrite("var_number", &var_columns_df::var_number)
        .def_readwrite("chrom", &var_columns_df::chrom)
        .def_readwrite("pos", &var_columns_df::pos)
        .def_readwrite("id", &var_columns_df::id)
        .def_readwrite("ref", &var_columns_df::ref)
        .def_readwrite("qual", &var_columns_df::qual)
        .def_readwrite("filter", &var_columns_df::filter)
        .def_readwrite("in_float", &var_columns_df::in_float)
        .def_readwrite("in_flag", &var_columns_df::in_flag)
        .def_readwrite("in_string", &var_columns_df::in_string)
        .def_readwrite("in_int", &var_columns_df::in_int)
        .def_readwrite("info_map1", &var_columns_df::info_map1)
        .def("print", &var_columns_df::print);

    py::class_<vcf_parsed>(m, "vcf_parsed")
        .def(py::init<>())
        .def_readwrite("id", &vcf_parsed::id)
        .def_readwrite("header", &vcf_parsed::header)
        .def_readwrite("INFO", &vcf_parsed::INFO)
        .def_readwrite("FORMAT", &vcf_parsed::FORMAT)
        .def_readwrite("var_columns", &vcf_parsed::var_columns)
        .def_readwrite("alt_columns", &vcf_parsed::alt_columns)
        .def_readwrite("samp_columns", &vcf_parsed::samp_columns)
        .def_readwrite("alt_sample", &vcf_parsed::alt_sample)
        .def("run", &vcf_parsed::run)
        .def("print_header", &vcf_parsed::print_header);

}
