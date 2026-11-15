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
#include "hydrodex/sound_shell_spectrum.hpp"
#include "hydrodex/constants.hpp"

namespace py = pybind11;
using namespace hydrodex;

using ProfileSolveMode = profile_calculator::ProfileSolveMode;
using NucleationHistory = sound_shell_model::NucleationHistory;
using ExpansionHistory = sound_shell_model::ExpansionHistory;

void bind_sound_shell_spectrum(pybind11::module_& m) {

    auto sound_shell_spectrum = py::class_<SoundShellSpectrum>(
        m, 
        "SoundShellSpectrum",
        "Class to calculate the velocity and gravitational wave spectrum in the"
        " sound shell model."
    );

    sound_shell_spectrum.def(
        py::init<
            const FluidProfile&,
            double, 
            NucleationHistory,
            ExpansionHistory,
            double,
            double,
            double
        >(), 
        py::arg("profile"), 
        py::arg("H_R"), 
        py::arg("nhist")=NucleationHistory::EXPONENTIAL, 
        py::arg("exhist")=ExpansionHistory::RADIATION_DOMINATION,
        py::arg("cs")=constants::cs_bag,
        py::arg("omega")=1.0/3.0,
        py::arg("Omega_star")=1.0,
        "Constructor for the SoundShellSpectrum from a FluidProfile instance."
        "\n\n"
        "Parameters\n"
        "----------\n"
        "profile : FluidProfile\n"
        "    FluidProfile instance\n"
        "H_R : float\n"
        "    Mean bubble separation at collision times Hubble rate, H* R*\n"
        "nhist : NucleationHistory, optional (default:"
        " NucleationHistory.EXPONENTIAL)\n    Nucleation history\n"
        "exhist : ExpansionHistory, optional (default:"
        " ExpansionHistory.RADIATION_DOMINATION)\n    Expansion history\n"
        "cs : float, optional (default: cs_bag)\n"
        "    Speed of sound cs. If negative or larger than 1, it is set to c-"
        " of profile\n"
        "omega : float, optional (default: 1/3)\n"
        "    Constant background equation-of-state parameter (only used if"
        " exhist is ExpansionHistory.CONSTANT_EOS)\n"
        "Omega_star : float, optional (default: 1)\n"
        "    Source-fluid energy-density fraction at the reference time"
    );

    sound_shell_spectrum.def(
        py::init([](
            double xi_w, 
            double alpha_N, 
            double H_R, 
            NucleationHistory nhist,
            ExpansionHistory exhist,
            double cs,
            double omega,
            double Omega_star,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ){
            return std::make_unique<SoundShellSpectrum>(
                xi_w, alpha_N, H_R, nhist, exhist, cs, omega, Omega_star, cp, cm,
                recover_dynamic_default(
                    step_size, profile_calculator::profile_step_size
                ), 
                recover_dynamic_default(
                    solve_mode, profile_calculator::profile_solve_mode
                )
            );
        }), 
        py::arg("xi_w"), 
        py::arg("alpha_N"), 
        py::arg("H_R"), 
        py::arg("nhist") = NucleationHistory::EXPONENTIAL, 
        py::arg("exhist")=ExpansionHistory::RADIATION_DOMINATION,
        py::arg("cs")=constants::cs_bag,
        py::arg("omega")=1.0/3.0,
        py::arg("Omega_star")=1.0,
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size") = py::none(),
        py::arg("solve_mode") = py::none(),
        "Constructor for the SoundShellSpectrum from the wall velocity xi_w and"
        " strength parameter alpha_N.\n\n"
        "Parameters\n"
        "----------\n"
        "xi_w : float\n"
        "    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Strength parameter alpha_N (far in front of the wall)\n"
        "H_R : float\n"
        "    Mean bubble separation at collision times Hubble rate, H* R*\n"
        "nhist : NucleationHistory, optional (default:"
        " NucleationHistory.EXPONENTIAL)\n    Nucleation history\n"
        "exhist : ExpansionHistory, optional (default:"
        " ExpansionHistory.RADIATION_DOMINATION)\n    Expansion history\n"
        "cs : float, optional (default: cs_bag)\n"
        "    Speed of sound cs in the fluid during the sound wave phase\n"
        "omega : float, optional (default: 1/3)\n"
        "    Constant background equation-of-state parameter (only used if"
        " exhist is ExpansionHistory.CONSTANT_EOS)\n"
        "Omega_star : float, optional (default: 1)\n"
        "    Source-fluid energy-density fraction at the reference time\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the bubble wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the bubble wall\n"
        "step_size : float, optional (default:"
        " profile_calculator.profile_step_size)\n    Step size for the fluid"
        " velocity (or coordinate) when solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved"
    );

    sound_shell_spectrum.def(
        "set_chi",
        py::overload_cast<double, double, size_t>(&SoundShellSpectrum::set_chi),
        py::arg("chi_min"), py::arg("chi_max"), py::arg("N"),
        "Set the values of chi = k T at which A^2 is evaluated.\n\n"
        "The values of chi = k T are spaced logarithmically with N points"
        " between (the logarithms of) chi_min and chi_max.\n\n"
        "Parameters\n"
        "----------\n"
        "chi_min : float\n"
        "    Minimal value of chi = k T\n"
        "chi_max : float\n"
        "    Maximal value of chi = k T\n"
        "N : int\n"
        "    Number of points"
    );

    sound_shell_spectrum.def(
        "set_chi",
        [](SoundShellSpectrum& self, py::array_t<double> chi_vec) { 
            self.set_chi(np2vec(chi_vec)); 
        },
        py::arg("chi_vec"),
        "Set the values of chi = k T at which A^2 is evaluated.\n\n"
        "Parameters\n"
        "----------\n"
        "chi_vec : numpy.ndarray\n"
        "    Array with values of chi = k T"
    );
    
    sound_shell_spectrum.def(
        "set_qR",
        py::overload_cast<double, double, size_t>(&SoundShellSpectrum::set_qR),
        py::arg("qR_min"), py::arg("qR_max"), py::arg("N"),
        "Set the values of the dimensionless momentum q R* at which the kinetic"
        " spectrum is evaluated.\n\n"
        "The values of q R* are spaced logarithmically with N points between"
        " (the logarithms of) qR_min and qR_max.\n\n"
        "Parameters\n"
        "----------\n"
        "qR_min : float\n"
        "    Minimal value of q R*\n"
        "qR_max : float\n"
        "    Maximal value of q R*\n"
        "N : int\n"
        "    Number of points"
    );

    sound_shell_spectrum.def(
        "set_qR",
        [](SoundShellSpectrum& self, py::array_t<double> qR_vec) { 
            self.set_qR(np2vec(qR_vec)); 
        },
        py::arg("qR_vec"),
        "Set the values of the dimensionless momentum q R* at which the kinetic"
        " spectrum is evaluated.\n\n"
        "Parameters\n"
        "----------\n"
        "qR_vec : numpy.ndarray\n"
        "    Array with values of q R*\n"
        "\nNote\n----\n"
        "The values must be logarithmically evenly spaced, since the kinetic"
        " spectrum is interpolated in log(q R*) when calculating the"
        " gravitational wave spectrum."
    );
    
    sound_shell_spectrum.def(
        "set_qR_from_kR",
        &SoundShellSpectrum::set_qR_from_kR,
        py::arg("N"),
        "Set the values of the dimensionless momentum q R* at which the kinetic"
        " spectrum is evaluated based on the values of k R* at which the"
        " gravitational wave spectrum is calculated.\n\n"
        "This function sets q R* making sure that the entire q R* range"
        " required to calculate the gravitational wave spectrum is covered."
        " It sets q R* to N logarithmically spaced points in the range"
        " [R* q(low), R* q(high)] with\n\n"
        "    q(low/high) = k(low/high) (1 -/+ cs)/(2 cs),\n\n"
        "where R* k(low/high) lie one logarithmic step below/above the lowest/"
        "highest values of k R*, with the step size determined by the spacing"
        " between the first two and last two values, respectively.\n\n"      
        "Parameters\n"
        "----------\n"
        "N : int\n"
        "    Number of points"
    );

    sound_shell_spectrum.def(
        "set_kR",
        [](
            SoundShellSpectrum& self, 
            double kR_min, 
            double kR_max, 
            size_t N, 
            py::object N_qR
        ) { self.set_kR(kR_min, kR_max, N, obj2idx(N_qR)); },
        py::arg("kR_min"),
        py::arg("kR_max"),
        py::arg("N"), 
        py::arg("N_qR")=py::none(),
        "Set the values of the dimensionless momentum k R* at which the"
        " gravitational wave spectrum is evaluated.\n\n"
        "The values of k R* are spaced logarithmically with N points between"
        " (the logarithms of) kR_min and kR_max.\n"
        "If N_qR is provided, the values for q R* for the evaluation of the"
        " velocity spectrum are set accordingly with N_qR points using"
        " set_qR_from_kR.\n\n"
        "Parameters\n"
        "----------\n"
        "kR_min : float\n"
        "    Minimal value of k R*\n"
        "kR_max : float\n"
        "    Maximal value of k R*\n"
        "N : int\n"
        "    Number of points (for k R*)\n"
        "N_qR : int or None, optional (default: None)\n"
        "    If provided, set_qR_from_kR is called to set q R* with N_qR points"
    );

    sound_shell_spectrum.def(
        "set_kR",
        [](
            SoundShellSpectrum& self, 
            py::array_t<double> kR_vec, 
            py::object N_qR
        ) { self.set_kR(np2vec(kR_vec), obj2idx(N_qR)); },
        py::arg("kR_vec"), py::arg("N_qR")=py::none(),
        "Set the values of the dimensionless momentum k R* at which the"
        " gravitational wave spectrum is evaluated.\n\n"
        "If N_qR is provided, the values for q R* for the evaluation of the"
        " velocity spectrum are set accordingly with N_qR points using"
        " set_qR_from_kR.\n\n"
        "Parameters\n"
        "----------\n"
        "kR_vec : numpy.ndarray\n"
        "    Array with values of k R*\n"
        "N_qR : int or None, optional (default: None)\n"
        "    If provided, set_qR_from_kR is called to set q R* with N_qR points"
    );

    sound_shell_spectrum.def(
        "calculate_A_squared",
        [](SoundShellSpectrum& self, double chi_min, double chi_max, size_t N) {
            auto [df,l] = self.calculate_A_squared(chi_min, chi_max, N);
            return py::make_tuple(vec2np(df), vec2np(l));
        },
        py::arg("chi_min"), py::arg("chi_max"), py::arg("N"),
        "Calculate the shape function A^2, Eq. (7.15) of [Barni+ (2027)].\n\n"
        "The shape function is evaluated on N logarithmically spaced values"
        " between chi_min and chi_max.\n"
        "The values of chi = k T are stored internally and can be retrieved"
        " calling get_chi.\n"
        "The values of A^2 are stored internally and can be retrieved calling"
        " get_shape_function.\n"
        "The function returns the vectors containing the values of the"
        " functions f'(chi) and l(chi) defined in Eqs. (7.16) of"
        " [Barni+ (2027)].\n\n"
        "Parameters\n"
        "----------\n"
        "chi_min : float\n"
        "    Minimal value of chi = k T\n"
        "chi_max : float\n"
        "    Maximal value of chi = k T\n"
        "N : int\n"
        "    Number of points\n\n"
        "Returns\n"
        "-------\n"
        "A2 : tuple(numpy.ndarray, numpy.ndarray)\n"
        "    Pair of vectors containing the values of f'(chi) and l(chi)"
    );

    sound_shell_spectrum.def(
        "calculate_A_squared",
        [](SoundShellSpectrum& self, py::array_t<double> chi_vec) {
            auto [df,l] = self.calculate_A_squared(np2vec(chi_vec));
            return py::make_tuple(vec2np(df), vec2np(l));
        },
        py::arg("chi_vec"), 
        "Calculate the shape function A^2, Eq. (7.15) of [Barni+ (2027)].\n\n"
        "The shape function is evaluated for the values provided in chi_vec.\n"
        "The values of chi = k T are stored internally and can be retrieved"
        " calling get_chi.\n"
        "The values of A^2 are stored internally and can be retrieved calling"
        " get_shape_function.\n"
        "The function returns the vectors containing the values of the"
        " functions f'(chi) and l(chi) defined in Eqs. (7.16) of"
        " [Barni+ (2027)].\n\n"
        "Parameters\n"
        "----------\n"
        "chi_vec : numpy.ndarray\n"
        "    Values of chi = k T at which the shape function is evaluated\n\n"
        "Returns\n"
        "-------\n"
        "A2 : tuple(numpy.ndarray, numpy.ndarray)\n"
        "    Pair of vectors containing the values of f'(chi) and l(chi)"
    );

    sound_shell_spectrum.def(
        "calculate_A_squared",
        [](SoundShellSpectrum& self) {
            auto [df,l] = self.calculate_A_squared();
            return py::make_tuple(vec2np(df), vec2np(l));
        },
        "Calculate the shape function A^2, Eq. (7.15) of [Barni+ (2027)].\n\n"
        "The shape function is evaluated for the values previously set using"
        " set_chi (or used in previous calls of calculate_A_squared) and can be"
        " retrieved calling get_chi.\n"
        "The values of A^2 are stored internally and can be retrieved calling"
        " get_shape_function.\n"
        "The function returns the vectors containing the values of the"
        " functions f'(chi) and l(chi) defined in Eqs. (7.16) of"
        " [Barni+ (2027)].\n\n"
        "Returns\n"
        "-------\n"
        "A2 : tuple(numpy.ndarray, numpy.ndarray)\n"
        "    Pair of vectors containing the values of f'(chi) and l(chi)"
    );

    sound_shell_spectrum.def(
        "calculate_kinetic_spectrum",
        py::overload_cast<double, double, size_t>(
            &SoundShellSpectrum::calculate_kinetic_spectrum
        ),
        py::arg("qR_min"), py::arg("qR_max"), py::arg("N"),
        "Calculate the kinetic spectrum (momentum dependence of the stationary"
        " velocity field UETC).\n\n"
        "The kinetic spectrum is evaluated on N logarithmically spaced values"
        " between qR_min and qR_max.\n"
        "The values of q R* are stored internally and can be retrieved calling"
        " get_qR.\n\n"
        "Parameters\n"
        "----------\n"
        "qR_min : float\n"
        "    Minimal value of q R*\n"
        "qR_max : float\n"
        "    Maximal value of q R*\n"
        "N : int\n"
        "    Number of points"
    );

    sound_shell_spectrum.def(
        "calculate_kinetic_spectrum",
        [](SoundShellSpectrum& self, py::array_t<double> qR_vec) { 
            self.calculate_kinetic_spectrum(np2vec(qR_vec)); 
        },
        py::arg("qR_vec"),
        "Calculate the kinetic spectrum (momentum dependence of the stationary"
        " velocity field UETC).\n\n"
        "The kinetic spectrum is evaluated for the values provided in qR_vec.\n"
        "The values of q R* are stored internally and can be retrieved calling"
        " get_qR.\n\n"
        "Parameters\n"
        "----------\n"
        "qR_vec : numpy.ndarray\n"
        "    Values of q R* at which the kinetic spectrum is evaluated\n"
        "\nNote\n----\n"
        "The values must be logarithmically evenly spaced, since the kinetic"
        " spectrum is interpolated in log(q R*) when calculating the"
        " gravitational wave spectrum."
    );
    
    sound_shell_spectrum.def(
        "calculate_kinetic_spectrum",
        py::overload_cast<>(&SoundShellSpectrum::calculate_kinetic_spectrum),
        "Calculate the kinetic spectrum (momentum dependence of the stationary"
        " velocity field UETC).\n\n"
        "The kinetic spectrum is evaluated for the values previously set using"
        " set_qR (or used in previous calls of calculate_kinetic_spectrum) and"
        " can be retrieved calling get_qR."
    );

    sound_shell_spectrum.def(
        "calculate_gw_spectrum",
        [](
            SoundShellSpectrum& self, 
            double kR_min, 
            double kR_max, 
            size_t N, 
            py::object N_qR
        ) { self.calculate_gw_spectrum(kR_min, kR_max, N, obj2idx(N_qR)); },
        py::arg("kR_min"), 
        py::arg("kR_max"), 
        py::arg("N"), 
        py::arg("N_qR")=py::none(),
        "Calculate the gravitational wave spectrum.\n\n"
        "The gravitational wave spectrum is evaluated on N logarithmically"
        " spaced values between kR_min and kR_max.\n"
        "The values of k R* are stored internally and can be retrieved calling"
        " get_kR.\n\n"
        "Parameters\n"
        "----------\n"
        "kR_min : float\n"
        "    Minimal value of k R*\n"
        "kR_max : float\n"
        "    Maximal value of k R*\n"
        "N : int\n"
        "    Number of points\n"
        "N_qR : int or None, optional (default: None)\n"
        "    If provided, set_qR_from_kR is called to set q R* with N_qR points"
    );

    sound_shell_spectrum.def(
        "calculate_gw_spectrum",
        [](
            SoundShellSpectrum& self, 
            py::array_t<double> kR_vec, 
            py::object N_qR
        ) { self.calculate_gw_spectrum(np2vec(kR_vec), obj2idx(N_qR)); },
        py::arg("kR_vec"), py::arg("N_qR")=py::none(),
        "Calculate the gravitational wave spectrum.\n\n"
        "The gravitational wave spectrum is evaluated for the values provided"
        " in kR_vec.\n"
        "The values of k R* are stored internally and can be retrieved calling"
        " get_kR.\n\n"
        "Parameters\n"
        "----------\n"
        "kR_vec : numpy.ndarray\n"
        "    Values of k R* at which the gravitational wave spectrum is"
        " evaluated\n"
        "N_qR : int or None, optional (default: None)\n"
        "    If provided, set_qR_from_kR is called to set q R* with N_qR points"
    );

    sound_shell_spectrum.def(
        "calculate_gw_spectrum",
        py::overload_cast<>(&SoundShellSpectrum::calculate_gw_spectrum),
        "Calculate the gravitational wave spectrum.\n\n"
        "The gravitational wave spectrum is evaluated for the values previously"
        " set using set_kR (or used in previous calls of calculate_gw_spectrum)"
        " and can be retrieved calling get_kR."
    );

    sound_shell_spectrum.def(
        "calculate_mean_squared_fluid_velocity", 
        &SoundShellSpectrum::calculate_mean_squared_fluid_velocity, 
        "Calculate the mean-square fluid velocity Uf^2 = <u^2>, cf. Eq. (7.17)"
        " of [Barni+ (2027)]."
    );

    sound_shell_spectrum.def(
        "calculate_anisotropic_stress_UETC_normalization", 
        &SoundShellSpectrum::calculate_anisotropic_stress_UETC_normalization, 
        "Calculate the normalization C of the anisotropic stress UETC."
    );
    
    sound_shell_spectrum.def(
        "calculate_anisotropic_stress_autocorrelator", 
        &SoundShellSpectrum::calculate_anisotropic_stress_autocorrelator, 
        "Calculate the anisotropic stress autocorrelator (momentum dependence"
        " of the anisotropic stress UETC).\n\n"
        "The anisotropic stress is evaluated at the same momenta as the"
        " gravitational wave spectrum, which can be set and retrieved calling"
        " set_kR and get_kR, respectively."
    );

    sound_shell_spectrum.def(
        "calculate_weighted_average_duration_dependence", 
        &SoundShellSpectrum::calculate_weighted_average_duration_dependence, 
        "Calculate the spectral-function-weighted average of the GW amplitude"
        " dependence on the source duration."
    );
    
    sound_shell_spectrum.def(
        "get_fluid_profile",
        &SoundShellSpectrum::get_fluid_profile,
        "Get the FluidProfile.\n\n"
        "Returns\n"
        "-------\n"
        "profile : FluidProfile\n"
        "    Fluid profile for which the gravitational wave spectrum is"
        " calculated"
    );
    
    sound_shell_spectrum.def(
        "get_chi",
        [](const SoundShellSpectrum& self) { return vec2np(self.get_chi()); },
        "Get the chi = k T vector at which the shape function A^2 is"
        " evaluated.\n\n"
        "Returns\n"
        "-------\n"
        "chi : numpy.ndarray\n"
        "    Values of chi = k T at which the shape function is evaluated"
    );
    
    sound_shell_spectrum.def(
        "get_qR",
        [](const SoundShellSpectrum& self) { return vec2np(self.get_qR()); },
        "Get the q R* vector at which the kinetic spectrum is evaluated.\n\n"
        "Returns\n"
        "-------\n"
        "qR : numpy.ndarray\n"
        "    Values of q R* at which the kinetic spectrum is evaluated"
    );
    
    sound_shell_spectrum.def(
        "get_kR",
        [](const SoundShellSpectrum& self) { return vec2np(self.get_kR()); },
        "Get the k R* vector at which the gravitational wave spectrum is"
        " evaluated.\n\n"
        "Returns\n"
        "-------\n"
        "kR : numpy.ndarray\n"
        "    Values of k R* at which the gravitational wave spectrum is"
        " evaluated"
    );
    
    sound_shell_spectrum.def(
        "get_shape_function",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_shape_function()); 
        },
        "Get the shape function A^2, Eq. (7.15) of [Barni+ (2027)].\n\n"
        "The corresponding values of chi can be obtained with get_chi.\n\n"
        "Returns\n"
        "-------\n"
        "A2 : numpy.ndarray\n"
        "    Shape function A^2"
    );
    
    sound_shell_spectrum.def(
        "get_kinetic_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_kinetic_spectrum()); 
        },
        "Get the dimensionless kinetic spectrum E_kin(q R*)/R* (i.e. the"
        " time-independent component of the velocity field UETC spectrum), cf."
        " Eq. (7.15) of [Barni+ (2027)].\n\n"
        "The corresponding values of q R* can be obtained with get_qR.\n\n"
        "Returns\n"
        "-------\n"
        "E_kin_over_R : numpy.ndarray\n"
        "    Dimensionless kinetic spectrum E_kin(q R*)/R*"
    );
    
    sound_shell_spectrum.def(
        "get_normalized_kinetic_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_normalized_kinetic_spectrum()); 
        },
        "Get the normalized kinetic spectrum zeta_kin(q R*) ="
        " E_kin(q)/(R* Uf^2), defined in Eq. (7.17) of [Barni+ (2027)].\n\n"
        "The corresponding values of q R* can be obtained with get_qR.\n\n"
        "Returns\n"
        "-------\n"
        "zeta_kin : numpy.ndarray\n"
        "    Normalized kinetic spectrum zeta_kin(q R*)"
    );
    
    sound_shell_spectrum.def(
        "get_plane_wave_velocity_spectral_density",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_plane_wave_velocity_spectral_density()); 
        },
        "Get the dimensionless spectral density of the plane wave components of"
        " the velocity field, S_v(q)/R*^3 = (2 pi^2 E_kin)/(q^2 R*^3), cf. Eq."
        " (4.17) of [Hindmarsh+ (2019)].\n\n"
        "The corresponding values of q R* can be obtained with get_qR.\n\n"
        "Returns\n"
        "-------\n"
        "S_v_over_R3 : numpy.ndarray\n"
        "    Dimensionless plane wave velocity spectral density S_v(q)/R*^3"
    );
    
    sound_shell_spectrum.def(
        "get_plane_wave_velocity_power_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_plane_wave_velocity_power_spectrum()); 
        },
        "Get the power spectrum of the plane wave components of the velocity"
        " field, P_v(q) = q^3 S_v(q) / (2 pi^2) = q E_kin(q).\n\n"
        "The corresponding values of q R* can be obtained with get_qR.\n\n"
        "Returns\n"
        "-------\n"
        "P_v : numpy.ndarray\n"
        "    Plane wave velocity power spectrum P_v(q)"
    );
    
    sound_shell_spectrum.def(
        "get_velocity_spectral_density",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_velocity_spectral_density()); 
        },
        "Get the dimensionless spectral density of the velocity field,"
        " S_v~(q)/R*^3 = 2 S_v(q)/R*^3 = (4 pi^2 E_kin)/(q^2 R*^3).\n\n"
        "The corresponding values of q R* can be obtained with get_qR.\n\n"
        "Returns\n"
        "-------\n"
        "S_v_tilde_over_R3 : numpy.ndarray\n"
        "    Dimensionless velocity field spectral density S_v~(q)/R*^3"
    );
    
    sound_shell_spectrum.def(
        "get_velocity_power_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_velocity_power_spectrum()); 
        },
        "Get the velocity power spectrum P_v~(q) = q^3 S_v~(q) / (2 pi^2) = 2 q"
        " E_kin(q), cf. Eq. (4.18) of [Hindmarsh+ (2019)].\n\n"
        "The corresponding values of q R* can be obtained with get_qR.\n\n"
        "Returns\n"
        "-------\n"
        "P_v_tilde : numpy.ndarray\n"
        "    Velocity field power spectrum P_v~(q)"
    );
    
    sound_shell_spectrum.def(
        "get_gravitational_wave_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_gravitational_wave_spectrum()); 
        },
        "Get the gravitational wave spectrum dOmega_GW*/dlog k at production,"
        " cf. Eq. (7.18) of [Barni+ (2027)].\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "Omega_GW : numpy.ndarray\n" 
        "    Gravitational wave spectrum Omega_GW(k)"
    );
    
    sound_shell_spectrum.def(
        "get_gravitational_wave_power_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_gravitational_wave_power_spectrum()); 
        },
        "Get the gravitational wave power spectrum P_gw(k), cf. Eq. (3.6) of"
        " [Hindmarsh+ (2019)].\n\n"
        "This is the same as get_gravitational_wave_spectrum.\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "P_gw : numpy.ndarray\n" 
        "    Gravitational wave power spectrum P_gw(k) = Omega_GW(k)"
    );
    
    sound_shell_spectrum.def(
        "get_gravitational_wave_spectral_density",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_gravitational_wave_spectral_density()); 
        },
        "Get the dimensionless gravitational wave spectral density"
        " S_gw(k)/R*^3.\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "S_gw_over_R3 : numpy.ndarray\n" 
        "    Dimensionless gravitational wave spectral density S_gw(k)/R*^3 ="
        " 2 pi^2 P_gw(k) / (k R*)^3"
    );
    
    sound_shell_spectrum.def(
        "get_anisotropic_stress_power_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_anisotropic_stress_power_spectrum()); 
        },
        "Get the anisotropic stress autocorrelator power spectrum P_Pi(k) ="
        " k E_Pi(k), cf. Eq. (7.10) of [Barni+ (2027)].\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "P_Pi : numpy.ndarray\n" 
        "    Anisotropic stress autocorrelator power spectrum P_Pi(k)"
    );
    
    sound_shell_spectrum.def(
        "get_normalized_anisotropic_stress_power_spectrum",
        [](SoundShellSpectrum& self) { return vec2np(
            self.get_normalized_anisotropic_stress_power_spectrum()
        ); },
        "Get the normalized anisotropic stress autocorrelator power spectrum"
        " (k R*)^3 zeta_Pi(k), cf. Eq. (7.25) of [Barni+ (2027)].\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "K3_zeta_Pi : numpy.ndarray\n" 
        "    Normalized anisotropic stress autocorrelator power spectrum"
        " (k R*)^3 zeta_Pi(k)"
    );
    
    sound_shell_spectrum.def(
        "get_anisotropic_stress_spectral_density",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_anisotropic_stress_spectral_density()); 
        },
        "Get the dimensionless anisotropic stress autocorrelator spectral"
        " density S_Pi(k)/R*^3 = 2 pi^2 E_Pi(k) / (k^2 R*^3).\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "S_Pi_over_R3 : numpy.ndarray\n" 
        "    Dimensionless anisotropic stress autocorrelator spectral density"
        " S_Pi(k)/R*^3"
    );
    
    sound_shell_spectrum.def(
        "get_normalized_gravitational_wave_spectrum",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_normalized_gravitational_wave_spectrum()); 
        },
        "Get the normalized gravitational wave spectrum (k R*)^3 zeta_GW"
        "(k R*), cf. Eq. (93) of [Roper+ (2024)].\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "K3_zeta_GW : numpy.ndarray\n" 
        "    Normalized gravitational wave spectrum (k R*)^3 zeta_GW(k R*)"
    );
    
    sound_shell_spectrum.def(
        "get_spectral_modification_function",
        [](SoundShellSpectrum& self) { 
            return vec2np(self.get_spectral_modification_function()); 
        },
        "Get the spectral modification function Delta~ = zeta_GW/zeta_Pi, cf."
        " Eq. (95) of [Roper+ (2024)].\n\n"
        "The corresponding values of k R* can be obtained with get_kR.\n\n"
        "Returns\n"
        "-------\n"
        "Delta_tilde : numpy.ndarray\n" 
        "    Spectral modification function Delta~"
    );

    sound_shell_spectrum.def(
        "get_mean_squared_fluid_velocity",
        &SoundShellSpectrum::get_mean_squared_fluid_velocity,
        "Get the mean-square fluid velocity Uf^2 = <u^2>, cf. Eq. (7.17) of"
        " [Barni+ (2027)].\n\n"
        "Returns\n"
        "-------\n"
        "Uf2 : float\n"
        "    Mean-square fluid velocity Uf^2"
    );
    
    sound_shell_spectrum.def(
        "get_kinetic_energy_density_parameter",
        &SoundShellSpectrum::get_kinetic_energy_density_parameter,
        "Get the total kinetic energy density parameter Omega_K ="
        " Omega* Gamma Uf^2, where Gamma is the adiabatic index of the fluid"
        " and Omega* the source-fluid energy-density fraction.\n\n"
        "Returns\n"
        "-------\n"
        "Omega_K : float\n"
        "    Kinetic energy density parameter Omega_K"
    );
    
    sound_shell_spectrum.def(
        "get_kinetic_spectrum_amplitude",
        &SoundShellSpectrum::get_kinetic_spectrum_amplitude,
        "Get the peak amplitude E_kin*/R* of the kinetic spectrum, cf. Eq. (41)"
        " of [Roper+ (2024)].\n\n"
        "Returns\n"
        "-------\n"
        "E_kin_peak_over_R : float\n"
        "    Peak amplitude of the kinetic spectrum in units of the bubble"
        " separation, E_kin*/R*"
    );
    
    sound_shell_spectrum.def(
        "get_kinetic_spectrum_peak_position",
        &SoundShellSpectrum::get_kinetic_spectrum_peak_position,
        "Get the peak position (q R*)_kin^peak of the kinetic spectrum.\n\n"
        "Returns\n"
        "-------\n"
        "qR_kin_peak : float\n"
        "    Peak position of the kinetic spectrum (q R*)_kin^peak"
    );

    sound_shell_spectrum.def(
        "get_normalized_kinetic_density_parameter",
        &SoundShellSpectrum::get_normalized_kinetic_density_parameter,
        "Get the normalized kinetic energy density parameter script_K ="
        " R* Uf^2 / (2 E_kin*), cf. Eqs. (42) and (43) of [Roper+ (2024)].\n\n"
        "Returns\n"
        "-------\n"
        "script_K : float\n"
        "    Kinetic density parameter normalized to the peak of the kinetic"
        " spectrum"
    );
    
    sound_shell_spectrum.def(
        "get_anisotropic_stress_UETC_normalization",
        &SoundShellSpectrum::get_anisotropic_stress_UETC_normalization,
        "Get the normalization script_C of the anisotropic stress spectral"
        " shape, cf. Eq. (7.25) of [Barni+ (2027)].\n\n"
        "Returns\n"
        "-------\n"
        "script_C : float\n"
        "    Normalization constant of the anisotropic stress spectral shape"
    );
    
    sound_shell_spectrum.def(
        "get_normalized_anisotropic_stress_UETC_amplitude",
        &SoundShellSpectrum::get_normalized_anisotropic_stress_UETC_amplitude,
        "Get the peak amplitude of the normalized anisotropic stress spectrum"
        " K^3 zeta_Pi(K).\n\n"
        "Returns\n"
        "-------\n"
        "K3_zeta_Pi_peak : float\n"
        "    Peak amplitude of the normalized anisotropic stress spectrum"
    );
    
    sound_shell_spectrum.def(
        "get_normalized_anisotropic_stress_UETC_peak_position",
        &SoundShellSpectrum::
            get_normalized_anisotropic_stress_UETC_peak_position,
        "Get the peak position (k R*)_Pi_peak of the normalized anisotropic"
        " stress spectrum.\n\n"
        "Returns\n"
        "-------\n"
        "kR_Pi_peak : float\n"
        "    Peak position of the normalized anisotropic stress spectrum"
    );
    
    sound_shell_spectrum.def(
        "get_weighted_average_duration_dependence",
        &SoundShellSpectrum::get_weighted_average_duration_dependence,
        "Get the spectral-function-weighted average Delta0~ of GW amplitude"
        " dependence on the source duration,"
        " Eq. (63) of [Roper+ (2024)].\n\n"
        "Returns\n"
        "-------\n"
        "Delta_0_tilde : float\n"
        "    Weighted average duration dependence Delta0~"
    );
    
    sound_shell_spectrum.def(
        "set_nucleation_history",
        &SoundShellSpectrum::set_nucleation_history,
        py::arg("nhist"),
        "Set the nucleation history type.\n\n"
        "Parameters\n"
        "----------\n"
        "nhist : NucleationHistory\n"
        "    Nucleation history (exponential or simultaneous)"
    );
    
    sound_shell_spectrum.def(
        "get_nucleation_history",
        &SoundShellSpectrum::get_nucleation_history,
        "Get the nucleation history type.\n\n"
        "Returns\n"
        "-------\n"
        "nhist : NucleationHistory\n"
        "    Nucleation history (exponential or simultaneous)"
    );
    
    sound_shell_spectrum.def(
        "recalculate_profile",
        [](
            SoundShellSpectrum& self, 
            double xi_w,
            double alpha_N,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ){
            return self.recalculate_profile(
                xi_w, alpha_N, cp, cm,
                recover_dynamic_default(
                    step_size, profile_calculator::profile_step_size
                ),
                recover_dynamic_default(
                    solve_mode, profile_calculator::profile_solve_mode
                )
            );
        },
        py::arg("xi_w"), 
        py::arg("alpha_N"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Change the wall velocity and/or strength parameter and recalculate"
        " the fluid profile.\n\n"
        "Parameters\n"
        "----------\n"
        "xi_w : float\n"
        "    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Strength parameter alpha_N (far in front of the wall)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the bubble wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the bubble wall\n"
        "step_size : float, optional (default:"
        " profile_calculator.profile_step_size)\n    Step size for the fluid"
        " velocity (or coordinate) when solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved"
    );
    
    sound_shell_spectrum.def(
        "get_wall_velocity",
        &SoundShellSpectrum::get_wall_velocity,
        "Get the value of the wall velocity xi_w.\n\n"
        "Returns\n"
        "-------\n"
        "xi_w : float\n"
        "    Wall velocity xi_w"
    );
    
    sound_shell_spectrum.def(
        "get_transition_strength",
        &SoundShellSpectrum::get_transition_strength,
        "Get the value of the strength parameter alpha_N.\n\n"
        "Returns\n"
        "-------\n"
        "alpha_N : float\n"
        "    Strength parameter alpha_N (far in front of the wall)"
    );
    
    sound_shell_spectrum.def(
        "set_expansion_history",
        &SoundShellSpectrum::set_expansion_history,
        py::arg("ehist"),
        "Set the expansion history type.\n\n"
        "Parameters\n"
        "----------\n"
        "ehist : ExpansionHistory\n"
        "    Expansion history\n"
        "\nRaises\n------\n"
        "ValueError\n    If ehist is ExpansionHistory.NONE"
    );

    sound_shell_spectrum.def(
        "get_expansion_history",
        &SoundShellSpectrum::get_expansion_history,
        "Get the expansion history type.\n\n"
        "Returns\n"
        "-------\n"
        "ehist : ExpansionHistory\n"
        "    Expansion history"
    );
    
    sound_shell_spectrum.def(
        "set_mean_bubble_separation",
        &SoundShellSpectrum::set_mean_bubble_separation,
        py::arg("H_R"),
        "Set the mean bubble separation H* R*.\n\n"
        "Parameters\n"
        "----------\n"
        "H_R : float\n"
        "    Mean bubble separation at collision times Hubble rate, H* R*"
    );
    
    sound_shell_spectrum.def(
        "get_mean_bubble_separation",
        &SoundShellSpectrum::get_mean_bubble_separation,
        "Get the mean bubble separation at collision times Hubble rate, H* R*."
        "\n\nReturns\n"
        "-------\n"
        "H_R : float\n"
        "    Mean bubble separation at collision times Hubble rate, H* R*"
    );
    
    sound_shell_spectrum.def(
        "set_source_duration",
        &SoundShellSpectrum::set_source_duration,
        py::arg("duration"),
        "Set the source duration to a fixed value.\n\n"
        "Parameters\n"
        "----------\n"
        "duration : float\n"
        "    Source duration H* delta_tau_fin in units of the Hubble time"
    );
    
    sound_shell_spectrum.def(
        "calculate_source_duration",
        &SoundShellSpectrum::calculate_source_duration,
        "Calculate the source duration.\n\n"
        "This function calculates the source duration from the mean bubble"
        " separation H* R* and the kinetic energy density Omega_K,"
        " H* delta_tau_fin = H* R* / sqrt(Omega_K / Omega*)."
    );
    
    sound_shell_spectrum.def(
        "get_source_duration",
        &SoundShellSpectrum::get_source_duration,
        "Get the duration of the sound wave source in units of the Hubble time."
        "\n\nReturns\n"
        "-------\n"
        "duration : float\n"
        "    Source duration H* delta_tau_fin in units of the Hubble time"
    );

    sound_shell_spectrum.def(
        "get_sound_speed",
        &SoundShellSpectrum::get_sound_speed,
        "Get the speed of sound in the fluid during the sound wave phase."
        "\n\nReturns\n"
        "-------\n"
        "cs : float\n    Speed of sound"
    );

    sound_shell_spectrum.def(
        "set_sound_speed",
        &SoundShellSpectrum::set_sound_speed,
        py::arg("cs"),
        "Set the speed of sound in the fluid during the sound wave phase.\n\n"
        "Parameters\n"
        "----------\n"
        "cs : float\n    Speed of sound"
    );

    sound_shell_spectrum.def(
        "get_omega",
        &SoundShellSpectrum::get_omega,
        "Get the constant background equation-of-state parameter.\n\n"
        "Returns\n"
        "-------\n"
        "omega : float\n    Constant background equation-of-state parameter"
    );

    sound_shell_spectrum.def(
        "set_omega",
        &SoundShellSpectrum::set_omega,
        py::arg("omega"),
        "Set the constant background equation-of-state parameter.\n\n"
        "Parameters\n"
        "----------\n"
        "omega : float\n    Constant background equation-of-state parameter"
    );

    sound_shell_spectrum.def(
        "get_Omega_star",
        &SoundShellSpectrum::get_Omega_star,
        "Get the source-fluid energy-density fraction at the reference time."
        "\n\nReturns\n"
        "-------\n"
        "Omega_star : float\n"
        "    Source-fluid energy-density fraction at the reference time"
    );

    sound_shell_spectrum.def(
        "set_Omega_star",
        &SoundShellSpectrum::set_Omega_star,
        py::arg("Omega_star"),
        "Set the source-fluid energy-density fraction at the reference time."
        "\n\nParameters\n"
        "----------\n"
        "Omega_star : float\n"
        "    Source-fluid energy-density fraction at the reference time"
    );
    
    sound_shell_spectrum.def(
        "set_correlator_integration_settings",
        &SoundShellSpectrum::set_correlator_integration_settings,
        py::arg("N_p"),
        py::arg("N_z"),
        py::arg("integration")=sound_shell_model::MomentumIntegration::FULL,
        "Set settings for the momentum integration of anisotropic stress"
        " correlators.\n\n"
        "Parameters\n"
        "----------\n"
        "N_p : int\n"
        "    Number of points in the integration over p R*\n"
        "N_z : int\n"
        "    Number of points in the integration over z\n"
        "integration : MomentumIntegration, optional (default:"
        " MomentumIntegration.FULL)\n"
        "    Approximation used in the momentum integration"
    );

    sound_shell_spectrum.def_static(
        "set_verbose",
        &SoundShellSpectrum::set_verbose,
        py::arg("b"),
        "Set verbosity.\n\n"
        "Parameters\n"
        "----------\n"
        "b : bool\n"
        "    Boolean value to which the verbosity flag is set"
    );

    sound_shell_spectrum.def_static(
        "is_verbose",
        &SoundShellSpectrum::is_verbose,
        "Get verbosity.\n\n"
        "Returns\n"
        "-------\n"
        "b : bool\n"
        "    Boolean value to which the verbosity flag is set"
    );
}