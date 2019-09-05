#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "Containers.h"
#include "SearchTypes/ConformationSearch.h"

namespace PNAB {
    std::string run(RuntimeParameters runtime_params, Backbone py_backbone,
                    std::vector<Base> py_bases, HelicalParameters hp, std::string prefix) {
        Backbone backbone(py_backbone.file_path, py_backbone.interconnects, py_backbone.linker, py_backbone.fixed_bonds);
        Bases bases(py_bases);

        ConformationSearch search(runtime_params, backbone, hp,  bases, prefix);
        std::string output = search.run();

        return output;

    }
}

namespace py = pybind11;

PYBIND11_MODULE(bind, m) {
    m.doc() = "Nucleic Acid Builder";

    py::class_<PNAB::RuntimeParameters>(m, "RuntimeParameters")
        .def(py::init())
        .def_readwrite("energy_filter", &PNAB::RuntimeParameters::energy_filter)
        .def_readwrite("max_distance", &PNAB::RuntimeParameters::max_distance)
        .def_readwrite("type", &PNAB::RuntimeParameters::type)
        .def_readwrite("search_algorithm", &PNAB::RuntimeParameters::search_algorithm)
        .def_readwrite("num_steps", &PNAB::RuntimeParameters::num_steps)
        .def_readwrite("dihedral_step", &PNAB::RuntimeParameters::dihedral_step)
        .def_readwrite("weighting_temperature", &PNAB::RuntimeParameters::weighting_temperature)
        .def_readwrite("monte_carlo_temperature", &PNAB::RuntimeParameters::monte_carlo_temperature)
        .def_readwrite("population_size", &PNAB::RuntimeParameters::population_size)
        .def_readwrite("mutation_rate", &PNAB::RuntimeParameters::mutation_rate)
        .def_readwrite("crossover_rate", &PNAB::RuntimeParameters::crossover_rate)
        .def_readwrite("strand", &PNAB::RuntimeParameters::strand)
        .def_readwrite("is_double_stranded", &PNAB::RuntimeParameters::is_double_stranded)
        .def_readwrite("is_hexad", &PNAB::RuntimeParameters::is_hexad)
        .def_readwrite("strand_orientation", &PNAB::RuntimeParameters::strand_orientation)
        ;

    py::class_<PNAB::HelicalParameters>(m, "HelicalParameters")
        .def(py::init())
        .def_readwrite("h_twist", &PNAB::HelicalParameters::h_twist)
        .def_readwrite("h_rise", &PNAB::HelicalParameters::h_rise)
        .def_readwrite("inclination", &PNAB::HelicalParameters::inclination)
        .def_readwrite("tip", &PNAB::HelicalParameters::tip)
        .def_readwrite("x_displacement", &PNAB::HelicalParameters::x_displacement)
        .def_readwrite("y_displacement", &PNAB::HelicalParameters::y_displacement)
        ;

    py::class_<PNAB::Backbone>(m, "Backbone")
        .def(py::init())
        .def_readwrite("file_path", &PNAB::Backbone::file_path)
        .def_readwrite("interconnects", &PNAB::Backbone::interconnects)
        .def_readwrite("linker", &PNAB::Backbone::linker)
        .def_readwrite("fixed_bonds", &PNAB::Backbone::fixed_bonds)
        ;        

    py::class_<PNAB::Base>(m, "Base")
        .def(py::init())
        .def_readwrite("name", &PNAB::Base::name)
        .def_readwrite("code", &PNAB::Base::code)
        .def_readwrite("pair_name", &PNAB::Base::pair_name)
        .def_readwrite("file_path", &PNAB::Base::file_path)
        .def_readwrite("linker", &PNAB::Base::linker)
        ;        

    m.def("run", &PNAB::run, "run");
}
