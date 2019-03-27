#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "Containers.h"
#include "SearchTypes/MonteCarloRotorSearch.h"

namespace PNAB {
    int run(RuntimeParameters runtime_params, Backbone py_backbone, std::vector<Base> py_bases, PyHelicalParameters py_hp) {
        runtime_params.validate();
        Backbone backbone(py_backbone.file_path, py_backbone.interconnects, py_backbone.linker);
        Bases bases(py_bases);
        HelicalParameters hp(py_hp);

        MonteCarloRotorSearch mcrs(runtime_params, backbone, hp,  bases);
        mcrs.run();

        return 0;

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
        .def_readwrite("parameter_file", &PNAB::RuntimeParameters::parameter_file)
        .def_readwrite("base_to_backbone_bond_length", &PNAB::RuntimeParameters::base_to_backbone_bond_length)
        .def_readwrite("num_steps", &PNAB::RuntimeParameters::num_steps)
        .def_readwrite("dihedral_discretization", &PNAB::RuntimeParameters::dihedral_discretization)
        .def_readwrite("angleStepSize", &PNAB::RuntimeParameters::angleStepSize)
        .def_readwrite("chain_length", &PNAB::RuntimeParameters::chain_length)
        .def_readwrite("algorithm", &PNAB::RuntimeParameters::algorithm)
        .def_readwrite("strand", &PNAB::RuntimeParameters::strand)
        .def_readwrite("is_double_stranded", &PNAB::RuntimeParameters::is_double_stranded)
        ;

    py::class_<PNAB::PyHelicalParameters>(m, "HelicalParameters")
        .def(py::init())
        .def_readwrite("tilt", &PNAB::PyHelicalParameters::tilt)
        .def_readwrite("roll", &PNAB::PyHelicalParameters::roll)
        .def_readwrite("twist", &PNAB::PyHelicalParameters::twist)
        .def_readwrite("shift", &PNAB::PyHelicalParameters::shift)
        .def_readwrite("slide", &PNAB::PyHelicalParameters::slide)
        .def_readwrite("rise", &PNAB::PyHelicalParameters::rise)
        .def_readwrite("buckle", &PNAB::PyHelicalParameters::buckle)
        .def_readwrite("propeller", &PNAB::PyHelicalParameters::propeller)
        .def_readwrite("opening", &PNAB::PyHelicalParameters::opening)
        .def_readwrite("shear", &PNAB::PyHelicalParameters::shear)
        .def_readwrite("stretch", &PNAB::PyHelicalParameters::stretch)
        .def_readwrite("stagger", &PNAB::PyHelicalParameters::stagger)
        .def_readwrite("inclination", &PNAB::PyHelicalParameters::inclination)
        .def_readwrite("tip", &PNAB::PyHelicalParameters::tip)
        .def_readwrite("x_displacement", &PNAB::PyHelicalParameters::x_displacement)
        .def_readwrite("y_displacement", &PNAB::PyHelicalParameters::y_displacement)
        ;

    py::class_<PNAB::Backbone>(m, "Backbone")
        .def(py::init())
        .def_readwrite("file_path", &PNAB::Backbone::file_path)
        .def_readwrite("interconnects", &PNAB::Backbone::interconnects)
        .def_readwrite("linker", &PNAB::Backbone::linker)
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
