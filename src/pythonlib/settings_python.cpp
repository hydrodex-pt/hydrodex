/*
 * Copyright (C) 2026 Giulio Barni, Eric Madge
 * This file is part of HydróDex.
 *
 * HydróDex is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * HydróDex is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * HydróDex. If not, see <https://www.gnu.org/licenses/>. 
 */

#include <pybind11/pybind11.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>

#include "hydrodex_bindings.hpp"
#include "python_helpers.hpp"
#include "hydrodex/settings.hpp"
#include "hydrodex/constants.hpp"

namespace py = pybind11;
using namespace hydrodex;

void bind_settings(pybind11::module_& m) {

    // set/get debug flag for additional command line prints
    m.def("get_debug", []() { return settings::debug; }, "Get the debug flag.");
    m.def(
        "set_debug",
        [](bool b) { settings::debug = b; },
        py::arg("b"),
        "Set the debug flag."
    );

    // set/get tolerance for considering two xi values as equal (jumps)
    m.def(
        "get_jump_tolerance",
        []() { return settings::jump_tolerance; },
        "Get the tolerance for considering two values of the coordinate as "
        "equal."
    );
    m.def(
        "set_jump_tolerance",
        [](double tol) {
            if ( !(tol >= 0.0) ) {
                throw std::invalid_argument(
                    "jump_tolerance must be non-negative, got "
                    + std::to_string(tol) + "."
                );
            }
            settings::jump_tolerance = tol;
        },
        py::arg("tol"),
        "Set the tolerance for considering two values of the coordinate as "
        "equal."
    );

    // set/get small velocity corresponding to zero when calculating profiles
    m.def(
        "get_v_zero",
        []() { return settings::v_zero; },
        "Get the small number corresponding to zero velocity when calculating "
        "profiles."
    );
    m.def(
        "set_v_zero",
        [](double v) {
            if ( !(v >= 0.0) ) {
                throw std::invalid_argument(
                    "v_zero must be non-negative, got " + std::to_string(v)
                    + "."
                );
            }
            settings::v_zero = v;
        },
        py::arg("v"),
        "Set the small number corresponding to zero velocity when calculating "
        "profiles."
    );

    py::native_enum<settings::CheckMode>(
        m, 
        "CheckMode", 
        "enum.Enum", 
        "Report mode options for non-critical consistency checks."
    )
        .value(
            "ERROR", 
            settings::CheckMode::ERROR, 
            "throw exception if a check fails"
        )
        .value(
            "WARNING", 
            settings::CheckMode::WARNING, 
            "show a warning if a check fails"
        )
        .value(
            "SILENT", 
            settings::CheckMode::SILENT, 
            "ignore failed checks"
        )
    .finalize();
    m.def(
        "get_consistency_check", 
        []() { return settings::consistency_check; }, 
        "Get the report mode for non-critical consistency checks."
    );
    m.def(
        "set_consistency_check", 
        [](settings::CheckMode mode) { settings::consistency_check = mode; }, 
        py::arg("mode"), 
        "Set the report mode for non-critical consistency checks."
    );

    m.def(
        "get_noninjectivity_check",
        []() { return settings::noninjectivity_check; },
        "Get the report mode for non-injectivity checks."
    );
    m.def(
        "set_noninjectivity_check",
        [](settings::CheckMode mode) {
            settings::noninjectivity_check = mode;
        },
        py::arg("mode"),
        "Set the report mode for non-injectivity checks."
    );

    // speed of sound
    m.attr("cs_bag")  = py::float_(constants::cs_bag);

}