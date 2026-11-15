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
#include <pybind11/numpy.h>

#include "hydrodex_bindings.hpp"
#include "python_helpers.hpp"
#include "hydrodex/fluid_profile.hpp"
#include "hydrodex/constants.hpp"

namespace py = pybind11;
using namespace hydrodex;

void bind_fluid_profile(pybind11::module_& m) {

    auto fluid_profile = py::class_<FluidProfile>(
        m, 
        "FluidProfile", 
        "Class representing the fluid profile of a bubble in a cosmological"
        " phase transition."
    );

    fluid_profile.def(
        py::init([](py::array_t<double> xi,
                    py::array_t<double> v,
                    py::array_t<double> w,
                    py::array_t<double> e,
                    py::array_t<double> p,
                    size_t idx_w,
                    py::object idx_sh,
                    double cp,
                    double cm
        ) {
            return std::make_unique<FluidProfile>( 
                np2vec(xi), 
                np2vec(v), 
                np2vec(w), 
                np2vec(e), 
                np2vec(p), 
                idx_w, 
                obj2idx(idx_sh),
                cp,
                cm
            );
        }),
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("w"), 
        py::arg("e"), 
        py::arg("p"),
        py::arg("idx_w"), 
        py::arg("idx_sh") = py::none(),
        py::arg("cp") = constants::cs_bag,
        py::arg("cm") = constants::cs_bag,
        "Constructor for the FluidProfile from precomputed profiles.\n\n"
        "Parameters\n"
        "----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Self-similar coordinate xi = r/t\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Fluid velocity profile v(xi)\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Enthalpy density profile w(xi)\n"
        "e : numpy.ndarray, shape (N,)\n"
        "    Energy density profile e(xi)\n"
        "p : numpy.ndarray, shape (N,)\n"
        "    Pressure profile p(xi)\n"
        "idx_w : int\n"
        "    Index of the first element outside the bubble\n"
        "idx_sh : int or None, optional (default: None)\n"
        "    Index of the first element outside the shock front\n"
        "    (use None if there is no shock front)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall"
    );

    fluid_profile.def(
        py::init([](
            double xi_w,
            double alpha_N,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ) {
            return std::make_unique<FluidProfile>( 
                xi_w, alpha_N, cp, cm, recover_dynamic_default(
                    step_size, profile_calculator::profile_step_size
                ), recover_dynamic_default(
                    solve_mode, profile_calculator::profile_solve_mode
                )
            );
        }),
        py::arg("xi_w"), 
        py::arg("alpha_N"), 
        py::arg("cp") = constants::cs_bag,
        py::arg("cm") = constants::cs_bag,
        py::arg("step_size") = py::none(),
        py::arg("solve_mode") = py::none(),
        "Constructor for the FluidProfile from the wall velocity xi_w and"
        " the pseudotrace phase transition strength parameter alphaN.\n\n"
        "Parameters\n"
        "----------\n"
        "xi_w : float\n"
        "    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Pseudotrace strength parameter alphaN far in front of the wall\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved"
    );

    fluid_profile.def("get_xi_coordinate", 
        [](const FluidProfile& self) {
            return vec2np(self.get_xi_coordinate()); 
        },
        "Get the self-similar coordinate xi=r/t of the profile.\n\n"
        "Returns\n"
        "-------\n"
        "xi : numpy.ndarray\n"
        "    Array containing the coordinate"
    );

    fluid_profile.def("get_fluid_velocity", 
        [](const FluidProfile& self) {
            return vec2np(self.get_fluid_velocity()); 
        },
        "Get the fluid velocity profile v(xi).\n\n"
        "Returns\n"
        "-------\n"
        "v : numpy.ndarray\n"
        "    Array containing the velocity"
    );

    fluid_profile.def("get_enthalpy_density", 
        [](const FluidProfile& self) {
            return vec2np(self.get_enthalpy_density()); 
        },
        "Get the enthalpy density profile w(xi).\n\n"
        "Returns\n"
        "-------\n"
        "w : numpy.ndarray\n"
        "    Array containing the enthalpy"
    );

    fluid_profile.def("get_energy_density", 
        [](const FluidProfile& self) {
            return vec2np(self.get_energy_density()); 
        },
        "Get the energy density profile e(xi).\n\n"
        "Returns\n"
        "-------\n"
        "e : numpy.ndarray\n"
        "    Array containing the energy"
    );

    fluid_profile.def("get_pressure", 
        [](const FluidProfile& self) {
            return vec2np(self.get_pressure()); 
        },
        "Get the pressure profile p(xi).\n\n"
        "Returns\n"
        "-------\n"
        "p : numpy.ndarray\n"
        "    Array containing the pressure"
    );

    fluid_profile.def("get_temperature",
        [](const FluidProfile& self, double Psi_N, double T_N) {
            return vec2np(self.get_temperature(Psi_N, T_N)); 
        },
        py::arg("Psi_N"), py::arg("T_N")=1.0,
        "Get the temperature profile T(xi).\n\n"
        "Parameters\n"
        "----------\n"
        "Psi_N : float\n"
        "    Enthalpy ratio Psi_N across the wall at nucleation\n"
        "T_N : float, optional (default: 1.0)\n"
        "    Nucleation temperature\n\n"
        "Returns\n"
        "-------\n"
        "T : numpy.ndarray\n"
        "    Array containing the temperature profile"
    );

    fluid_profile.def("get_entropy_density",
        [](const FluidProfile& self, double Psi_N, double s_N) {
            return vec2np(self.get_entropy_density(Psi_N, s_N)); 
        },
        py::arg("Psi_N"), py::arg("s_N")=1.0,
        "Get the entropy density profile s(xi).\n\n"
        "Parameters\n"
        "----------\n"
        "Psi_N : float\n"
        "    Enthalpy ratio Psi_N across the wall at nucleation\n"
        "s_N : float, optional (default: 1.0)\n"
        "    Entropy density at nucleation\n\n"
        "Returns\n"
        "-------\n"
        "s : numpy.ndarray\n"
        "    Array containing the entropy density profile"
    );


    fluid_profile.def("get_wall_index", 
        &FluidProfile::get_wall_index,
        "Get the index of the element in xi corresponding to the bubble wall"
        " position.\n\n"
        "Returns\n"
        "-------\n"
        "idx : int\n"
        "    Index of the bubble wall position (first element outside the wall)"
    );

    fluid_profile.def("get_shock_index",
        [](const FluidProfile& self) -> pybind11::object { 
            return idx2obj(self.get_shock_index()); 
        },
        "Get the index of the element in xi corresponding to the shock front"
        " position.\n\n"
        "Returns\n"
        "-------\n"
        "idx : int or None\n    Index of the shock front position (first"
        " element outside the shock), or None if there is no shock front"
    );

    fluid_profile.def("get_wall_velocity",
        &FluidProfile::get_wall_velocity,
        "Get the wall velocity xi_w.\n\n"
        "Returns\n"
        "-------\n"
        "xi_w : float\n"
        "    Value of the wall velocity xi_w"
    );

    fluid_profile.def("get_transition_strength",
        &FluidProfile::get_transition_strength,
        "Get the pseudotrace phase transition strength alphaN.\n\n"
        "Returns\n"
        "-------\n"
        "alpha_N : float\n"
        "    Value of the transition strength alphaN"
    );

    fluid_profile.def("get_sound_speed_plus",
        &FluidProfile::get_sound_speed_plus,
        "Get the speed of sound c+ in the false vacuum phase.\n\n"
        "Returns\n"
        "-------\n"
        "cp : float\n"
        "    Speed of sound c+ in the false vacuum phase"
    );

    fluid_profile.def("get_sound_speed_minus",
        &FluidProfile::get_sound_speed_minus,
        "Get the speed of sound c- in the true vacuum phase.\n\n"
        "Returns\n"
        "-------\n"
        "cm : float\n"
        "    Speed of sound c- in the true vacuum phase"
    );

    fluid_profile.def("get_v_plus",
        &FluidProfile::get_v_plus,
        "Get the wall-frame fluid velocity v+ ahead of the wall.\n\n"
        "Returns\n"
        "-------\n"
        "v_plus : float\n"
        "    Wall-frame fluid velocity v+ ahead of the wall"
    );

    fluid_profile.def("get_v_minus",
        &FluidProfile::get_v_minus,
        "Get the wall-frame fluid velocity v- behind the wall.\n\n"
        "Returns\n"
        "-------\n"
        "v_minus : float\n"
        "    Wall-frame fluid velocity v- behind the wall"
    );

    fluid_profile.def("get_wall_transition_strength",
        &FluidProfile::get_wall_transition_strength,
        "Get the pseudotrace phase transition strength alpha+ ahead of the"
        " wall.\n\n"
        "Returns\n"
        "-------\n"
        "alpha_plus : float\n"
        "    Pseudotrace phase transition strength alpha+ ahead of the bubble"
        " wall"
    );

    fluid_profile.def("get_shock_position",
        [](const FluidProfile& self) -> pybind11::object {
            if ( !self.has_shock() ) { return py::none(); }
            return py::float_(self.get_shock_position());
        },
        "Get the self-similar coordinate xi_sh of the shock front.\n\n"
        "Returns\n"
        "-------\n"
        "xi_sh : float or None\n"
        "    Position xi_sh of the shock front, or None if there is no shock"
        " front"
    );

    fluid_profile.def("get_solution_type",
        &FluidProfile::get_solution_type,
        "Get the type of the hydrodynamic solution.\n\n"
        "Returns\n"
        "-------\n"
        "sol_type : HydroSolutionType\n"
        "    Hydrodynamic solution type"
    );

    fluid_profile.def("has_shock",
        &FluidProfile::has_shock,
        "Check whether the profile has a shock front.\n\n"
        "Returns\n"
        "-------\n"
        "b : bool\n"
        "    True if the profile has a shock front, False otherwise"
    );

    fluid_profile.def("get_mean_enthalpy_density",
        &FluidProfile::get_mean_enthalpy_density,
        "Get the mean enthalpy density (far in front of the wall).\n\n"
        "Returns\n"
        "-------\n"
        "w_bar : float\n"
        "    Mean enthalpy density w_bar"
    );

    fluid_profile.def(
        "get_efficiency",
        &FluidProfile::get_efficiency,
        py::arg("normalization") =
            profile_calculator::EfficiencyNormalization::TotalEnergyDensity,
        "Get the efficiency K (or kappa_sw).\n\n"
        "Parameters\n"
        "----------\n"
        "normalization : EfficiencyNormalization, optional (default:"
        " profile_calculator.EfficiencyNormalization.TotalEnergyDensity)\n"
        "    EfficiencyNormalization switch to determine how to normalize the"
        " efficiency\n\n"
        "Returns\n"
        "-------\n"
        "kappa_or_K : float\n"
        "    Efficiency K or kappa_sw"
    );

    fluid_profile.def_static("set_verbose",
        &FluidProfile::set_verbose,
        py::arg("b"),
        "Set verbosity.\n\n"
        "Parameters\n"
        "----------\n"
        "b : bool\n"
        "    Boolean value to which the verbosity flag is set"
    );

    fluid_profile.def_static("is_verbose",
        &FluidProfile::is_verbose,
        "Get verbosity.\n\n"
        "Returns\n"
        "-------\n"
        "b : bool\n"
        "    Boolean value to which the verbosity flag is set"
    );
}