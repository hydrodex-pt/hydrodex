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

#include "hydrodex_bindings.hpp"
#include "python_helpers.hpp"
#include "hydrodex/sound_shell_model.hpp"
#include "hydrodex/constants.hpp"

namespace py = pybind11;
using namespace hydrodex;
using namespace hydrodex::sound_shell_model;
using MomentumIntegration 
    = AnisotropicStressCorrelatorSettings::MomentumIntegration;

void bind_sound_shell_model(pybind11::module_& m) {

    // Nucleation history types
    py::native_enum<NucleationHistory>(
        m, "NucleationHistory", "enum.Enum", "Nucleation history options"
    )
        .value(
            "EXPONENTIAL",  
            NucleationHistory::EXPONENTIAL,  
            "exponential nucleation"
        )
        .value(
            "SIMULTANEOUS", 
            NucleationHistory::SIMULTANEOUS, 
            "simultaneous nucleation"
        )
    .finalize();

    // Expansion history types
    py::native_enum<ExpansionHistory>(
        m, "ExpansionHistory", "enum.Enum", "Expansion history options"
    )
        .value(
            "CONSTANT_EOS",
            ExpansionHistory::CONSTANT_EOS,
            "assume a constant equation-of-state parameter"
        )
        .value(
            "RADIATION_DOMINATION", 
            ExpansionHistory::RADIATION_DOMINATION, 
            "assume a radiation dominated universe"
        )
        .value(
            "FLAT_UNIVERSE", 
            ExpansionHistory::FLAT_UNIVERSE,
            "assume a flat (non-expanding) universe"
        )
        .value(
            "NONE", 
            ExpansionHistory::NONE, 
            "no expansion history / time integration (i.e. use Delta=1 to"
            " calculate the autocorrelation zeta_Pi)"
        )
    .finalize();

    // Approximation used in the momentum integration
    py::native_enum<MomentumIntegration>(
        m, 
        "MomentumIntegration", 
        "enum.Enum", 
        "Approximation used in the momentum integration of the anisotropic"
        " stress correlator"
    )
        .value(
            "FULL",
            MomentumIntegration::FULL,
            "Full momentum integration, using Eq. (7.21) of [Barni+ (2027)]"
            " for Delta."
        )
        .value(
            "LOW_FREQUENCY",
            MomentumIntegration::LOW_FREQUENCY,
            "Use the low-frequency limit Delta_0 or Delta_0_flat, i.e. Eq."
            " (7.23) of [Barni+ (2027)] for Delta."
        )
        .value(
            "INFINITE_DURATION",
            MomentumIntegration::INFINITE_DURATION,
            "Use the infinite-duration (high-frequency) limit of Delta,"
            " k * delta_tau_fin >> 1, i.e. Eq. (7.24) of [Barni+ (2027)]. In"
            " this limit, Delta reduces to Dirac delta-functions, simplifying"
            " the momentum integration."
        )
    .finalize();

    py::class_<AnisotropicStressCorrelatorSettings> (
        m, 
        "AnisotropicStressCorrelatorSettings",
        "Settings for the momentum and time integration of the anisotropic"
        " stress and gravitational wave spectrum.\n"
        "\nNote\n----\nThe parameter omega is ignored unless universe is set"
        " to ExpansionHistory.CONSTANT_EOS."
    )
        .def(
            py::init<>(),
            "Default constructor (for calculating the anisotropic stress"
            " autocorrelation power spectrum k E_Pi(k))."
        )
        .def(
            py::init<size_t, size_t>(),
            py::arg("N_p"), 
            py::arg("N_z"),
            "Constructor for settings for calculating the anisotropic stress"
            " autocorrelation power spectrum k E_Pi(k).\n\n"
            "Note: This sets universe to ExpansionHistory.NONE, freq_approx to"
            " MomentumIntegration.FULL.\n"
            "\nParameters\n----------\n"
            "N_p : int\n    Number of points in the integration over"
            " P = p R*\n"
            "N_z : int\n    Number of points in the integration over z"
        )
        .def(
            py::init<
                double, 
                double, 
                ExpansionHistory, 
                MomentumIntegration, 
                double,
                double,
                double
            >(),
            py::arg("H_delta_tau"), 
            py::arg("H_R"),
            py::arg("universe")=ExpansionHistory::RADIATION_DOMINATION,
            py::arg("freq_approx")=MomentumIntegration::FULL,
            py::arg("c_s") = constants::cs_bag,
            py::arg("omega") = 1.0/3.0,
            py::arg("Omega_star") = 1.0,
            "Constructor for time integration settings for the gravitational"
            " wave spectrum Omega_GW.\n"
            "\nParameters\n----------\n"
            "H_delta_tau : float\n    Source duration in units of the Hubble"
            " time, H* delta_tau_fin\n"
            "H_R : float\n    Mean bubble separation in units of the Hubble"
            " size, H* R*\n"
            "universe : ExpansionHistory, optional (default:"
            " ExpansionHistory.RADIATION_DOMINATION)\n"
            "    Assumptions used in the time integration\n"
            "freq_approx : MomentumIntegration, optional (default:"
            " MomentumIntegration.FULL)\n"
            "    Approximation used in the momentum integration of the time"
            " integrals\n"
            "c_s : float, optional (default: cs_bag)\n"
            "    Speed of sound cs in the stable phase\n"
            "omega : float, optional (default: 1.0/3.0)\n"
            "    Equation-of-state parameter omega of the Universe\n"
            "Omega_star : float, optional (default: 1.0)\n"
            "    Fractional energy density Omega* of the fluid (kinetic+"
            "thermal) at production"
        )
        .def(
            py::init<
                size_t, 
                size_t, 
                double, 
                double, 
                ExpansionHistory,
                MomentumIntegration,
                double,
                double,
                double
            >(),
            py::arg("N_p"), 
            py::arg("N_z"),
            py::arg("H_delta_tau"), 
            py::arg("H_R"),
            py::arg("universe")=ExpansionHistory::RADIATION_DOMINATION,
            py::arg("freq_approx")=MomentumIntegration::FULL,
            py::arg("c_s") = constants::cs_bag,
            py::arg("omega") = 1.0/3.0,
            py::arg("Omega_star") = 1.0,
            "Full constructor.\n"
            "\nParameters\n----------\n"
            "N_p : int\n"
            "    Number of points in the integration over P = p R*\n"
            "N_z : int\n"
            "    Number of points in the integration over z\n"
            "H_delta_tau : float\n    Source duration in units of the Hubble"
            " time, H* delta_tau_fin\n"
            "H_R : float\n"
            "    Mean bubble separation in units of the Hubble size, H* R*\n"
            "universe : ExpansionHistory, optional (default: "
            "ExpansionHistory.RADIATION_DOMINATION)\n"
            "    Assumptions used in the time integration\n"
            "freq_approx : MomentumIntegration, optional (default:"
            " MomentumIntegration.FULL)\n    Approximation used in the momentum"
            " integration of the time integrals\n"
            "c_s : float, optional (default: cs_bag)\n"
            "    Speed of sound cs in the stable phase\n"
            "omega : float, optional (default: 1.0/3.0)\n"
            "    Equation-of-state parameter omega of the Universe\n"
            "Omega_star : float, optional (default: 1.0)\n"
            "    Fractional energy density Omega* of the fluid (kinetic+"
            "thermal) at production"
        )
        .def_readwrite(
            "N_p", 
            &AnisotropicStressCorrelatorSettings::N_p, 
            "Number of points in the integration over p R*"
            " (default: sound_shell_model.N_pR_default)"
        )
        .def_readwrite(
            "N_z", 
            &AnisotropicStressCorrelatorSettings::N_z, 
            "Number of points in the integration over z"
            " (default: sound_shell_model.N_z_default)"
        )
        .def_readwrite(
            "H_delta_tau", 
            &AnisotropicStressCorrelatorSettings::H_delta_tau, 
            "Source duration in units of the Hubble time, H* delta_tau_fin"
        )
        .def_readwrite(
            "H_R", 
            &AnisotropicStressCorrelatorSettings::H_R, 
            "Mean bubble separation in units of the Hubble size, H* R*"
        )
        .def_readwrite(
            "universe", 
            &AnisotropicStressCorrelatorSettings::universe, 
            "Assumptions used in the time integration"
        )
        .def_readwrite(
            "freq_approx", 
            &AnisotropicStressCorrelatorSettings::freq_approx,
            "Approximation used in the momentum integration of the time"
            " integrals"
        )
        .def_readwrite(
            "c_s", 
            &AnisotropicStressCorrelatorSettings::c_s, 
            "Speed of sound cs in the stable phase"
        )
        .def_readwrite(
            "omega", 
            &AnisotropicStressCorrelatorSettings::omega, 
            "Equation-of-state parameter of the Universe"
        )
        .def_readwrite(
            "Omega_star",
            &AnisotropicStressCorrelatorSettings::Omega_star,
            "Fractional energy density Omega* of the fluid (kinetic + thermal)"
            " at production"
        )
    ;  

    m.def(
        "adiabatic_index",
        &adiabatic_index,
        py::arg("cs") = constants::cs_bag,
        "Get the mean adiabatic index Gamma = w/e.\n\n"
        "This calculates the mean adiabatic index from the speed of sound,"
        " Gamma = 1 + cs^2.\n"
        "\nParameters\n----------\n"
        "cs : float, optional (default: cs_bag)\n"
        "    Speed of sound cs\n"
        "\nReturns\n"
        "-------\n"
        "Gamma : float\n"
        "    Mean adiabatic index Gamma = w/e"
    );

    m.def(
        "hubble_times_tau",
        &hubble_times_tau,
        py::arg("omega") = 1.0/3.0,
        "Comoving Hubble rate times conformal time H* tau*.\n\n"
        "This function returns the relation between the comoving Hubble rate"
        " and the conformal time for a given equation-of-state parameter omega"
        " of the Universe, H tau = 2/(1 + 3 omega).\n"
        "\nParameters\n----------\n"
        "omega : float, optional (default: 1/3)\n"
        "    Equation-of-state parameter omega of the Universe\n"
        "\nReturns\n-------\n"
        "H_tau : float\n    Comoving Hubble rate times conformal time H tau"
    );

    m.def(
        "nu_exp",
        &nu_exp,
        py::arg("T"),
        "Bubble lifetime distribution function for exponential nucleation.\n\n"
        "This is Eq. (4.27) of [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "T : float\n    Dimensionless bubble lifetime T~ = beta T_i\n"
        "\nReturns\n-------\n"
        "nu_exp : float\n    Bubble lifetime distribution evaluated at T for"
        " exponential nucleation"
    );

    m.def(
        "nu_sim",
        &nu_sim,
        py::arg("T"),
        "Bubble lifetime distribution function for simultaneous nucleation.\n\n"
        "This is Eq. (4.32) of [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "T : float\n    Dimensionless bubble lifetime T~ = beta T_i\n"
        "\nReturns\n-------\n"
        "nu_sim : float\n    Bubble lifetime distribution evaluated at T for"
        " simultaneous nucleation"
    );

    m.def(
        "calculate_sine_transforms",
        [](
            py::array_t<double> chi_arr,
            py::array_t<double> xi,
            py::array_t<double> v,
            py::array_t<double> w,
            py::array_t<double> e,
            size_t idx_wall,
            py::object idx_sh,
            double c_s
        ) {
            std::vector<double> A2, df, l;
            std::tie(A2, df, l) = calculate_sine_transforms(
                np2vec(chi_arr), 
                np2vec(xi), 
                np2vec(v), 
                np2vec(w), 
                np2vec(e), 
                idx_wall, 
                obj2idx(idx_sh),
                c_s
            );
            return py::make_tuple(vec2np(A2), vec2np(df), vec2np(l));
        },
        py::arg("chi_arr"),
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("w"), 
        py::arg("e"), 
        py::arg("idx_wall"), 
        py::arg("idx_sh")=py::none(),
        py::arg("c_s")=constants::cs_bag,
        "Calculate the (derivative of) sine transforms of the velocity and"
        " energy spectrum.\n\n"
        "This function calculates the derivative of the sine transform f(chi)"
        " of the velocity spectrum v, the sine transform l(chi) of the energy"
        " fluctuation variable lambda = (e-e_bar)/w_bar, and the absolute value"
        " squared of their sum,"
        " ``|A(chi)|^2 = 1/4 [(f'(chi))^2 + (cs l(chi))^2]``,"
        " given in Eqs. (7.15) and (7.16) of [Barni+ (2027)].\n"
        "\nParameters\n----------\n"
        "chi_arr : numpy.ndarray, shape (N,)\n"
        "    Array of values of the dimensionless wavenumber chi = k Tn\n"
        "xi : numpy.ndarray, shape (M,)\n"
        "    Values of the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (M,)\n    Fluid velocity profile v(xi)\n"
        "w : numpy.ndarray, shape (M,)\n    Enthalpy density profile w(xi)\n"
        "e : numpy.ndarray, shape (M,)\n    Energy density profile e(xi)\n"
        "idx_wall : int\n    Index of the bubble wall position"
        " (first element outside the wall)\n"
        "idx_sh : int or None, optional (default: None)\n"
        "    Index of the shock front position (first element outside the"
        " shock, use None if there is no shock front)\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound cs in the stable phase\n"
        "\nReturns\n-------\n"
        "A2 : numpy.ndarray, shape (N,)\n    Shape function ``|A(chi)|^2``\n"
        "df : numpy.ndarray, shape (N,)\n    Derivative f'(chi) of the sine"
        " transform of the velocity spectrum\n"
        "l : numpy.ndarray, shape (N,)\n"
        "    Sine transform l(chi) of the energy fluctuation variable"
    );

    m.def(
        "scaled_kinetic_spectrum",
        [](
            py::array_t<double> kR, 
            py::array_t<double> chi_arr, 
            py::array_t<double> A2_arr, 
            double xi_wall, 
            NucleationHistory nuc_hist
        ) {
            return vec2np(scaled_kinetic_spectrum(
                np2vec(kR), np2vec(chi_arr), np2vec(A2_arr), xi_wall, nuc_hist
            ));
        },
        py::arg("kR"), 
        py::arg("chi_arr"), 
        py::arg("A2_arr"), 
        py::arg("xi_wall"), 
        py::arg("nuc_hist")=NucleationHistory::EXPONENTIAL,
        "Scaled kinetic spectrum Ekin(k)/R*.\n\n"
        "This is the momentum dependent coefficient Ekin(k) of the velocity"
        " field unequal-time correlator,\n\n"
        "    Ekin(k, tau1, tau2) = Ekin(k) cos(k cs (tau1-tau2))\n\n"
        "cf. Eq. (7.15) of [Barni+ (2027)]. We use the dimensionless momentum"
        " z = k R* as input variable and divide by R* to make the kinetic"
        " spectrum dimensionless.\n\n"
        "Note that this is related to the velocity spectral density of the"
        " plane wave components of the velocity field, Eq. (4.17) of"
        " [Hindmarsh+ (2019)], P_v(k) = 2 pi^2 Ekin(k)/k^2.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectrum is evaluated\n"
        "chi_arr : numpy.ndarray, shape (M,)\n"
        "    Values of chi = T q = T~ q/beta at which ``|A|^2`` is provided\n"
        "A2_arr : numpy.ndarray, shape (M,)\n    Precomputed values of the"
        " shape function ``|A|^2`` used in the integration\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "nuc_hist : NucleationHistory, optional (default:"
        " NucleationHistory.EXPONENTIAL)\n    Nucleation history (exponential"
        " or simultaneous) used for the collision time distribution\n"
        "\nReturns\n-------\n"
        "Ekin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k R*)/R*"
    );

    m.def(
        "kinetic_spectral_shape",
        [](py::array_t<double> E_kin_over_R, double Uf2) { 
            return vec2np(kinetic_spectral_shape(np2vec(E_kin_over_R), Uf2)); 
        },
        py::arg("E_kin_over_R"),
        py::arg("Uf2"),
        "Spectral shape zeta_kin(k R*) of the kinetic spectrum.\n\n"
        "This is the spectral shape zeta_kin(k R*) of the kinetic spectrum,"
        " Ekin(k) = Uf^2 R* zeta_kin(k R*), where Uf is the root-mean-square"
        " fluid velocity and R* is the mean bubble separation. Cf. Eq. (7.17)"
        " of [Barni+ (2027)].\n"
        "\nNote\n----\n"
        "Our definition here differs from Eq. (41) of [Roper+ (2024)], where"
        " the spectral shape is instead normalized to the maximum amplitude"
        " Ekin* of the kinetic spectrum.\n"
        "\nParameters\n----------\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k)/R*\n"
        "Uf2 : float\n"
        "    Mean-square fluid velocity Uf^2\n"
        "\nReturns\n-------\n"
        "zeta_kin : numpy.ndarray, shape (N,)\n"
        "    Kinetic spectral shape zeta_kin(k R*) = Ekin(k)/R*/Uf^2"
    );

    m.def(
        "mean_squared_fluid_velocity",
        [](py::array_t<double> kR, py::array_t<double> E_kin_over_R) { 
            return mean_squared_fluid_velocity(
                np2vec(kR), np2vec(E_kin_over_R)
            ); 
        },
        py::arg("kR"), py::arg("E_kin_over_R"),
        "Mean-square fluid velocity Uf^2.\n\n"
        "This is the root-mean-square fluid velocity squared, cf. Eq. (7.17)"
        " of [Barni+ (2027)], calculated from the kinetic spectrum (or"
        " velocity-field UETC at equal times, see also Eqs. (4.33) and (4.34)"
        " of [Hindmarsh+ (2019)]).\n"
        "\nNote\n----\nThe sound shell model assumes non-relativistic fluid"
        " velocities and a constant entropy density, w(x) = w_bar ~ wN, where"
        " wN is the enthalpy density at the nucleation temperature. Thus, the"
        " root-mean-square fluid velocity coincides with the"
        " enthalpy-weighted average fluid velocity.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectrum is evaluated\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n    Scaled kinetic spectrum"
        " Ekin(k)/R*\n"
        "\nReturns\n-------\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2"
    );

    m.def(
        "mean_squared_fluid_velocity",
        [](
            py::array_t<double> chi_arr, 
            py::array_t<double> A2_arr, 
            double xi_wall
        ) { 
            return mean_squared_fluid_velocity(
                np2vec(chi_arr), np2vec(A2_arr), xi_wall
            ); 
        },
        py::arg("chi_arr"), py::arg("A2_arr"), py::arg("xi_wall"),
        "Mean-square fluid velocity Uf^2.\n\n"
        "This is the root-mean-square fluid velocity squared, cf. Eqs. (4.33)"
        " and (4.34) of [Hindmarsh+ (2019)].\n"
        "The integral over T~ in Eq. (4.33) has been carried out analytically,"
        " assuming that the nucleation history is either exponential or"
        " simultaneous [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "chi_arr : numpy.ndarray, shape (N,)\n"
        "    Values of chi = T q = T~ q/beta at which ``|A|^2`` is provided\n"
        "A2_arr : numpy.ndarray, shape (N,)\n    Precomputed values of the"
        " shape function ``|A|^2`` used in the integration\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "\nReturns\n-------\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2"
    );

    m.def(
        "kinetic_energy_density_parameter",
        py::overload_cast<double,double,double>(
            &kinetic_energy_density_parameter
        ),
        py::arg("Uf2"), py::arg("Gamma")=4.0/3.0, py::arg("Omega_star")=1.0,
        "Total kinetic energy density parameter OmegaK.\n\n"
        "This is the total kinetic energy density OmegaK = rho_kin/rho_crit as"
        " a fraction of the critical energy density, where the kinetic energy"
        " is given by rho_kin = <w gamma^2 v^2>. It is related to the"
        " root-mean-square fluid velocity by OmegaK = Omega* Gamma Uf^2, where"
        " Gamma = w_bar/e_bar is the mean adiabatic index of the fluid in the"
        " stable phase, and Omega* = rho_*/rho_crit is the total fractional"
        " energy density of the fluid (kinetic + thermal). The latter is set to"
        " Omega* = 1 if the fluid dominates the Universe. Cf. Eq. (7.17) of"
        " [Barni+ (2027)].\n"
        "\nNote\n----\n"
        "This assumes Gamma=4/3 by default. To obtain OmegaK = 1/2 Uf^2 as"
        " defined in Eq. (42) of [Roper+ (2024)], set Gamma to 0.5.\n"
        "\nParameters\n----------\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2\n"
        "Gamma : float, optional (default: 4.0/3.0)\n    Mean adiabatic index"
        " Gamma = w_bar/e_bar of the fluid in the stable phase\n"
        "Omega_star : float, optional (default: 1.0)\n    Total fractional"
        " energy density Omega* of the fluid (kinetic + thermal)\n"
        "\nReturns\n-------\n"
        "OmegaK : float\n    Kinetic energy density parameter OmegaK\n"
    );

    m.def(
        "kinetic_energy_density_parameter",
        [](
            py::array_t<double> kR, 
            py::array_t<double> E_kin_over_R, 
            double Gamma,
            double Omega_star
        ) { 
            return kinetic_energy_density_parameter(
                np2vec(kR), np2vec(E_kin_over_R), Gamma, Omega_star
            );
        },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("Gamma")=4.0/3.0,
        py::arg("Omega_star")=1.0,
        "Total kinetic energy density parameter OmegaK.\n\n"
        "This is the total kinetic energy density OmegaK = rho_kin/rho_crit as"
        " a fraction of the critical energy density, calculated from the"
        " kinetic spectrum (or velocity-field UETC at equal times, cf. Eq. (42)"
        " of [Roper+ (2024)]). Cf. Eq. (7.17) of [Barni+ (2027)].\n"
        "\nNote\n----\nSet Gamma to 0.5 and Omega_star to 1.0 to obtain"
        " OmegaK as defined in Eq. (42) of [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectrum is evaluated\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k)/R*\n"
        "Gamma : float, optional (default: 4.0/3.0)\n    Mean adiabatic index"
        " Gamma = w_bar/e_bar of the fluid in the stable phase\n"
        "Omega_star : float, optional (default: 1.0)\n    Total fractional"
        " energy density Omega* of the fluid (kinetic + thermal)\n"
        "\nReturns\n-------\n"
        "OmegaK : float\n    Kinetic energy density parameter OmegaK"
    );

    m.def(
        "kinetic_energy_density_parameter",
        [](
            double xi_wall,
            py::array_t<double> chi_arr, 
            py::array_t<double> A2_arr, 
            double Gamma,
            double Omega_star
        ) { 
            return kinetic_energy_density_parameter(
                xi_wall, np2vec(chi_arr), np2vec(A2_arr), Gamma, Omega_star
            );
        },
        py::arg("xi_wall"), 
        py::arg("chi_arr"), 
        py::arg("A2_arr"), 
        py::arg("Gamma")=4.0/3.0,
        py::arg("Omega_star")=1.0,
        "Total kinetic energy density parameter OmegaK.\n\n"
        "This is the total kinetic energy density OmegaK = rho_kin/rho_crit as"
        " a fraction of the critical energy density.\n"
        "The integral over T~ in Eq. (4.33) has been carried out analytically,"
        " assuming that the nucleation history is either exponential or"
        " simultaneous, as in the step from Eq. (4.33) to (4.34) in"
        " [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "chi_arr : numpy.ndarray, shape (N,)\n"
        "    Values of chi = T q = T~ q/beta at which ``|A|^2`` is provided\n"
        "A2_arr : numpy.ndarray, shape (N,)\n    Precomputed values of the"
        " shape function ``|A|^2`` used in the integration\n"
        "Gamma : float, optional (default: 4.0/3.0)\n    Mean adiabatic index"
        " Gamma = w_bar/e_bar of the fluid in the stable phase\n"
        "Omega_star : float, optional (default: 1.0)\n    Total fractional"
        " energy density Omega* of the fluid (kinetic + thermal)\n"
        "\nReturns\n-------\n"
        "OmegaK : float\n    Kinetic energy density parameter OmegaK"
    );

    m.def(
        "anisotropic_stress_power_spectrum_amplitude",
        [](py::array_t<double> kR, py::array_t<double> zeta_kin) { 
            return anisotropic_stress_power_spectrum_amplitude(
                np2vec(kR), np2vec(zeta_kin)
            ); 
        },
        py::arg("kR"), py::arg("zeta_kin"),
        "Normalized amplitude script_C of the anisotropic stress power spectrum"
        " k E_Pi.\n\nThis is the constant script_C defined in Eq. (7.25) of"
        " [Barni+ (2027)]. It describes the normalized amplitude of the"
        " anisotropic stress power spectrum [cf. also Eqs. (44) and (46) of"
        " [Roper+ (2024)]]\n"
        "\n    k E_Pi(tau_1,tau_2,k) = 2 K^3 Uf^4 script_C zeta_Pi(tau-, K)\n\n"
        "where K = k R* and Uf is the RMS fluid velocity.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectral shape is evaluated\n"
        "zeta_kin : numpy.ndarray, shape (N,)\n"
        "    Kinetic spectral shape zeta_kin(k R*)\n"
        "\nReturns\n-------\n"
        "script_C : float\n    Normalized amplitude script_C of the anisotropic"
        " stress power spectrum k E_Pi"
    );

    m.def(
        "anisotropic_stress_correlator_integration",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R, 
            AnisotropicStressCorrelatorSettings settings
        ) {
            return vec2np(anisotropic_stress_correlator_integration(
                np2vec(kR), E_kin_over_R, settings
            ));
        },
        py::arg("kR"), py::arg("E_kin_over_R"), py::arg("settings"),
        "Calculate the nested momentum integral of the anisotropic stress UETC"
        " E_Pi or gravitational wave spectrum OmegaGW*.\n\n"
        "This function calculates the nested momentum (and time) integrals for"
        " different choices of the time dependence function Delta. The choice"
        " of Delta and the normalization are controlled by the universe and"
        " freq_approx parameters of the AnisotropicStressCorrelatorSettings"
        " settings.\n"
        "If settings.universe is set to ExpansionHistory.NONE, the function"
        " returns the anisotropic stress (auto-correlation) power spectrum"
        " k E_Pi(k), i.e. Eq. (7.10) of [Barni+ (2027)].\n"
        "For ExpansionHistory.RADIATION_DOMINATION,"
        " ExpansionHistory.CONSTANT_EOS, or ExpansionHistory.FLAT_UNIVERSE, the"
        " function returns the gravitational wave spectrum OmegaGW(k) at"
        " production, i.e. Eq. (7.18) of [Barni+ (2027)] setting T_GW=1,"
        " using the approximation specified in settings.freq_approx.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the anisotropic stress UETC is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R*)\n"
        "settings : AnisotropicStressCorrelatorSettings\n"
        "    AnisotropicStressCorrelatorSettings object containing the settings"
        " for the momentum and time integration\n"
        "\nReturns\n-------\n"
        "K3_zeta : numpy.ndarray, shape (N,)\n"
        "    Integrated anisotropic stress UETC or GW spectrum"
    );

    m.def(
        "anisotropic_stress_power_spectrum",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R, 
            py::object N_p, 
            py::object N_z
        ) {
            return vec2np(anisotropic_stress_power_spectrum(
                np2vec(kR), 
                E_kin_over_R, 
                recover_dynamic_default(N_p, N_pR_default),
                recover_dynamic_default(N_z, N_z_default)
            ));
        },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("N_p")=py::none(), 
        py::arg("N_z")=py::none(),
        "Anisotropic stress auto-correlation power spectrum k E_Pi.\n\n"
        "This is the power spectrum of the anisotropic stress auto-correlation,"
        " k E_Pi(k) = k E_Pi(tau,tau,k), i.e. Eq. (7.10) of [Barni+ (2027)],"
        " corresponding to the integration in"
        " anisotropic_stress_correlator_integration"
        " with Delta = 1 (settings.universe = ExpansionHistory.NONE).\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the anisotropic stress UETC is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R*)\n"
        "N_p : int, optional (default: sound_shell_model.N_pR_default)\n"
        "    Number of points in the integration over p R*\n"
        "N_z : int, optional (default: sound_shell_model.N_z_default)\n"
        "    Number of points in the integration over z\n"
        "\nReturns\n-------\n"
        "k_E_Pi : numpy.ndarray, shape (N,)\n    Anisotropic stress"
        " auto-correlator power spectrum k E_Pi(tau,tau,k)"
    );

    m.def(
        "anisotropic_stress_power_spectral_shape",
        [](
            py::array_t<double> k_E_Pi, 
            double script_C, 
            double Uf2
        ) {
            return vec2np(anisotropic_stress_power_spectral_shape(
                np2vec(k_E_Pi), script_C, Uf2
            ));
        },
        py::arg("k_E_Pi"), 
        py::arg("script_C"), 
        py::arg("Uf2"), 
        "Spectral shape of the anisotropic stress UETC power spectrum (k R*)^3"
        " zeta_Pi(k R*).\n\n"
        "This is the combination K^3 zeta_Pi(K) with script_C zeta_Pi(K)"
        " defined in Eq. (7.25) of [Barni+ (2027)]. It describes the spectral"
        " shape of the power spectrum of the anisotropic stress UETC k E_Pi."
        " The constant script_C can be calculated with"
        " anisotropic_stress_power_spectrum_amplitude(qR, zeta_kin).\n"
        "\nParameters\n----------\n"
        "k_E_Pi : numpy.ndarray, shape (N,)\n"
        "    Anisotropic stress auto-correlator k E_Pi\n" 
        "script_C : float\n    Amplitude parameter script_C of the anisotropic"
        " stress power spectrum\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2\n"
        "\nReturns\n-------\n"
        "K3_zeta_Pi : numpy.ndarray, shape (N,)\n    Spectral shape of the"
        " anisotropic stress UETC power spectrum (k R*)^3 zeta_Pi(k R*)"
    );  

    m.def(
        "infinite_duration_gw_power_spectrum",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R, 
            double Upsilon_H_R, 
            double c_s,
            double Omega_star,
            py::object N_p
        ) {
            return vec2np(infinite_duration_gw_power_spectrum(
                np2vec(kR), E_kin_over_R, Upsilon_H_R, c_s, Omega_star, 
                recover_dynamic_default(N_p, N_pR_default)
            ));
        },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("Upsilon_H_R"), 
        py::arg("c_s")=constants::cs_bag,
        py::arg("Omega_star")=1.0,
        py::arg("N_p")=py::none(), 
        "Approximated gravitational wave power spectrum OmegaGW* in the"
        " infinite duration limit.\n\n"
        "This is the gravitational wave power spectrum OmegaGW* at production"
        " [cf. anisotropic_stress_correlator_integration] in the infinite"
        " source duration limit, Eq. (B3) of [Roper+ (2024)] (with T_GW=1),"
        " cf. also Eq. (3.13) of [Giombi+ (2025)]."
        " In the k >> 1/delta_tau_fin limit, the time integration factor Delta"
        " becomes a Dirac-delta distribution, cf. Eq. (7.24) of [Barni+"
        " (2027)], simplifying the momentum integration.\n\n"
        "In a flat universe, this corresponds to P_GW in Eq. (3.48) of"
        " [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the spectrum is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R*)\n"
        "Upsilon_H_R : float\n"
        "    Finite source lifetime suppression factor, Upsilon, times the mean"
        " bubble separation normalized to the Hubble radius, H* R*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound cs in the stable phase\n"
        "Omega_star : float, optional (default: 1)\n"
        "    Fractional energy density Omega* of the fluid"
        " (kinetic + thermal)\n"
        "N_p : int, optional (default: sound_shell_model.N_pR_default)\n"
        "    Number of points in the integration over p R*\n"
        "\nReturns\n-------\n"
        "OmegaGW : numpy.ndarray, shape (N,)\n    GW power spectrum OmegaGW* at"
        " production in the k delta_tau_fin >> 1 limit"
    );

    m.def(
        "gravitational_wave_power_spectrum",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R,
            double H_delta_tau,
            double H_R,
            ExpansionHistory expansion,
            MomentumIntegration integration,
            double c_s,
            double omega,
            double Omega_star,
            py::object N_p,
            py::object N_z
        ) { return vec2np(gravitational_wave_power_spectrum(
            np2vec(kR), 
            E_kin_over_R,
            H_delta_tau, 
            H_R, 
            expansion, 
            integration,  
            c_s,
            omega,
            Omega_star,
            recover_dynamic_default(N_p, N_pR_default), 
            recover_dynamic_default(N_z, N_z_default)
        )); },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("H_delta_tau"), 
        py::arg("H_R"), 
        py::arg("expansion")=ExpansionHistory::RADIATION_DOMINATION,
        py::arg("integration")=MomentumIntegration::FULL, 
        py::arg("c_s")=constants::cs_bag,
        py::arg("omega")=1.0/3.0,
        py::arg("Omega_star")=1.0,
        py::arg("N_p")=py::none(), 
        py::arg("N_z")=py::none(),
        "Gravitational wave power spectrum OmegaGW*(k) at production.\n\n"
        "This is Eq. (7.26) of [Barni+ (2027)].\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the spectrum is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R*)\n"
        "H_delta_tau : float\n"
        "    Source duration in units of the Hubble time, H* delta_tau_fin\n"
        "H_R : float\n"
        "    Mean bubble separation in units of the Hubble size, H* R*\n"
        "expansion : ExpansionHistory, optional (default:"
        " ExpansionHistory.RADIATION_DOMINATION)\n"
        "    Assumed ExpansionHistory\n"
        "integration : MomentumIntegration, optional (default:"
        " MomentumIntegration.FULL)\n    Method for calculating the momentum"
        " integrals in anisotropic_stress_correlator_integration\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound cs in the stable phase\n"
        "omega : float, optional (default: 1/3)\n"
        "    Equation-of-state parameter omega of the Universe\n"
        "Omega_star : float, optional (default: 1)\n"
        "    Fractional energy density Omega* of the fluid"
        " (kinetic + thermal)\n"
        "N_p : int, optional (default: sound_shell_model.N_pR_default)\n"
        "    Number of points in the integration over p R*\n"
        "N_z : int, optional (default: sound_shell_model.N_z_default)\n"
        "    Number of points in the integration over z\n"
        "\nReturns\n-------\n"
        "OmegaGW : numpy.ndarray, shape (N,)\n"
        "    GW power spectrum OmegaGW*(k) at production"
    );

    m.def(
        "time_integration_exponent",
        py::overload_cast<double,double>(&time_integration_exponent),
        py::arg("omega"), py::arg("c_s"),
        "Acoustic-velocity damping exponent in the time integration.\n\n"
        "This is the parameter sigma_bar = (1 - 3 cs^2)/(1 + 3 omega)"
        " describing the damping exponent of the acoustic fluid velocity, cf."
        " Eq. (7.3) of [Barni+ (2027)]. It is not the exponent of the final"
        " gravitational-wave time kernel.\n"
        "\nParameters\n----------\n"
        "omega : float\n"
        "    Equation-of-state parameter omega of the Universe\n"
        "c_s : float\n"
        "    Speed of sound cs in the stable phase\n"
        "\nReturns\n-------\n"
        "sigma_bar : float\n    Acoustic-velocity damping exponent sigma_bar"
    );

    m.def(
        "time_integration_exponent",
        py::overload_cast<double>(&time_integration_exponent),
        py::arg("omega"),
        "Exponent modifier in the time integration.\n\n"
        "This is the parameter sigma = (1 - 3 omega)/(1 + 3 omega) modifying"
        " the exponent 1-sigma in the time integrals due to the modified"
        " equation of state of the Universe, i.e. Eq. (7.3) of [Barni+"
        " (2027)].\n"
        "\nParameters\n----------\n"
        "omega : float\n"
        "    Equation-of-state parameter omega of the Universe\n"
        "\nReturns\n-------\n"
        "sigma : float\n    Exponent modifier sigma in the time integration"
    );

    m.def(
        "Delta_rad",
        &Delta_rad,
        py::arg("H_delta_tau"), 
        py::arg("k_tau"), 
        py::arg("p_tau"), 
        py::arg("pt_tau"),
        py::arg("c_s")=constants::cs_bag,
        "Time integral function Delta(delta_tau_fin, k, p, p~) for a radiation"
        " dominated universe.\n\n"
        "This is the function Delta that contains the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC (cf. Eq. (52) of [Roper+ (2024)]),"
        " assuming a radiation dominated universe and that the fluid red-shifts"
        " like radiation.\n\n"
        "We use the dimensionless variables q tau* = q/H* for the wave numbers"
        " q=k,p,p~ and H delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we assume radiation domination, we have H* = 1/tau*,"
        " so that q tau* = (q R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "k_tau : float\n    Dimensionless wave number k tau*\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "pt_tau : float\n    Dimensionless wave number p~ tau*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "\nReturns\n-------\n"
        "Delta : float\n"
        "    Value of the time integral Delta(delta_tau_fin, k, p, p~)"
    );

    m.def(
        "Delta_mn",
        &Delta_mn,
        py::arg("H_delta_tau"), py::arg("pmn_tau"),
        "Intermediate function Delta_mn(delta_tau_fin, p_mn) for the time"
        " integral in a radiation-dominated universe.\n\n"
        "This is the intermediate function Delta_mn(delta_tau_fin, p_mn) used"
        " in the calculation of the time integral Delta_rad, cf. Eq. (56) of"
        " [Roper+ (2024)].\n\n"
        "We use the dimensionless variables pmn tau* = pmn/H* and H"
        " delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "pmn_tau : float\n    Dimensionless wave number p_mn tau*\n"
        "\nReturns\n-------\n"
        "Delta_mn : float\n    Value of the intermediate time integral"
        " contribution Delta_mn(delta_tau_fin, p_mn)"
    );

    m.def(
        "Delta_flat",
        &Delta_flat,
        py::arg("H_delta_tau"),
        py::arg("k_tau"), 
        py::arg("p_tau"), 
        py::arg("pt_tau"),
        py::arg("c_s")=constants::cs_bag,
        "Time integral function Delta(delta_tau_fin, k, p, p~) for a flat"
        " universe.\n\n"
        "This is the function Delta that contains the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC (cf. Eq. (52) of [Roper+ (2024)]),"
        " assuming a flat universe.\n\n"
        "We use the dimensionless variables q tau* = q/H* for the wave numbers"
        " q=k,p,p~ and H delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we assume radiation domination, we have H* = 1/tau*,"
        " so that q tau* = (q R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "k_tau : float\n    Dimensionless wave number k tau*\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "pt_tau : float\n    Dimensionless wave number p~ tau*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "\nReturns\n-------\n"
        "Delta : float\n    Value of the time integral"
        " Delta_flat(delta_tau_fin, k, p, p~) for a flat universe"
    );

    m.def(
        "Delta_mn_flat",
        &Delta_mn_flat,
        py::arg("H_delta_tau"), py::arg("pmn_tau"),
        "Intermediate function Delta_mn_flat(delta_tau_fin, p_mn) for the time"
        " integral in a flat universe.\n\n"
        "This is Eq. (59) of [Roper+ (2024)].\n\n"
        "We use the dimensionless variables pmn tau* = pmn/H* and H"
        " delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "pmn_tau : float\n    Dimensionless wave number p_mn tau*\n"
        "\nReturns\n-------\n"
        "Delta_mn : float\n"
        "    Value of the time integral contribution Delta_mn_flat"
        "(delta_tau_fin, p_mn) neglecting the expansion of the Universe"
    );

    m.def(
        "Delta_const_eos",
        &Delta_const_eos,
        py::arg("H_delta_tau"),
        py::arg("k_over_H"),
        py::arg("p_over_H"),
        py::arg("pt_over_H"),
        py::arg("c_s")=constants::cs_bag,
        py::arg("omega")=1.0/3.0,
        "Time integral function Delta(delta_tau_fin, k, p, p~) for a universe"
        " with a constant equation-of-state parameter.\n\n"
        "This is the function Delta that contains the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC, Eqs. (7.19) and (7.21) of [Barni+"
        " (2027)], assuming a universe with a constant equation-of-state"
        " parameter as well as a fluid with a constant speed of sound.\n\n"
        "We use the dimensionless variables q/H* for the wave numbers q=k,p,p~"
        " and H* delta_tau_fin as inputs.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "k_over_H : float\n    Dimensionless wave number k/H*\n"
        "p_over_H : float\n    Dimensionless wave number p/H*\n"
        "pt_over_H : float\n    Dimensionless wave number p~/H*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "omega : float, optional (default: 1/3)\n"
        "    Equation-of-state parameter omega of the Universe\n"
        "\nReturns\n-------\n"
        "Delta : float\n"
        "    Value of the time integral Delta(delta_tau_fin, k, p, p~)"
    );

    m.def(
        "Delta_mn_const_eos",
        &Delta_mn_const_eos,
        py::arg("delta_tau_over_tau"), py::arg("pmn_tau"), py::arg("sigma"),
        "Intermediate function Delta_mn(delta_tau_fin, p_mn) for the time"
        " integral in a universe with a constant equation-of-state parameter.\n"
        "\nThis is the intermediate function Delta_mn(delta_tau_fin, p_mn) used"
        " in the calculation of the time integral Delta_const_eos, Eq. (7.21)"
        " of [Barni+ (2027)].\n\n"
        "We use the dimensionless variables pmn tau* and delta_tau_fin/tau* as"
        " inputs. This function assumes that the Universe has a constant"
        " equation of state and that the fluid has a constant speed of sound.\n"
        "\nParameters\n----------\n"
        "delta_tau_over_tau : float\n    Duration of the gravitational wave"
        " source normalized to the initial time delta_tau_fin/tau*\n"
        "pmn_tau : float\n    Dimensionless wave number p_mn tau*\n"
        "sigma : float\n    Exponent modifier sigma in the time integration\n"
        "\nReturns\n-------\n"
        "Delta_mn : float\n    Value of the intermediate time integral"
        " contribution Delta_mn"
    );

    m.def(
        "Delta_low_freq_rad",
        &Delta_low_freq_rad,
        py::arg("H_delta_tau"), 
        py::arg("p_tau"), 
        py::arg("c_s")=constants::cs_bag,
        "Time integral function in the low frequency limit"
        " Delta0(delta_tau_fin, p) = lim(k->0) Delta(delta_tau_fin,k,p,p~).\n\n"
        "This is the low-frequency limit Delta0 of the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC, cf. Eq. (61) of [Roper+ (2024)],"
        " assuming radiation domination.\n\n"
        "We use the dimensionless variables p tau* = p/H* and H delta_tau_fin"
        " = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we assume radiation domination, we have H* = 1/tau*,"
        " so that p tau* = (p R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "\nReturns\n-------\n"
        "Delta0 : float\n"
        "    Value of the time integral Delta0 in the low frequency limit"
    );

    m.def(
        "Delta_low_freq_flat",
        &Delta_low_freq_flat,
        py::arg("H_delta_tau"), 
        py::arg("p_tau"),
        py::arg("c_s")=constants::cs_bag,
        "Time integral function in the low frequency limit in a flat"
        " universe, Delta0_flat(delta_tau_fin, p) = lim(k->0)"
        " Delta_flat(delta_tau_fin,k,p,p~).\n\n"
        "This is the low-frequency limit Delta0_flat of the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC, cf. Eq. (61) of [Roper+ (2024)],"
        " assuming a flat universe.\n\n"
        "We use the dimensionless variables p tau* = p/H* and H delta_tau_fin"
        " = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we assume radiation domination, we have"
        " H* = 1/tau*, so that p tau* = (p R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "\nReturns\n-------\n"
        "Delta0 : float\n    Value of the time integral Delta0_flat in the low"
        " frequency limit in a flat universe"
    );

    m.def(
        "Delta_low_freq_const_eos",
        &Delta_low_freq_const_eos,
        py::arg("H_delta_tau"),
        py::arg("p_over_H"),
        py::arg("c_s") = constants::cs_bag,
        py::arg("omega") = 1.0/3.0,
        "Time integral function in the low frequency limit in a universe with"
        " a constant equation of state,"
        " Delta0(delta_tau_fin, p) = lim(k->0) Delta(delta_tau_fin,k,p,p~).\n\n"
        "This is the low-frequency limit Delta0 of the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC, Eq. (7.23) of [Barni+ (2027)],"
        " assuming a universe with a constant equation-of-state parameter and"
        " a constant speed of sound for the fluid.\n\n"
        "We use the dimensionless variables p/H* and H* delta_tau_fin as"
        " inputs.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "p_over_H : float\n    Dimensionless wave number p/H*\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "omega : float, optional (default: 1/3)\n"
        "    Equation-of-state parameter omega of the Universe\n"
        "\nReturns\n-------\n"
        "Delta0 : float\n"
        "    Value of the time integral Delta0 in the low frequency limit"
    );

    m.def(
        "weighted_average_Delta_0",
        [](
            py::array_t<double> kR, 
            py::array_t<double> E_kin_over_R,  
            double H_delta_tau, 
            double H_R, 
            ExpansionHistory expansion,
            double c_s,
            double omega
        ) { return weighted_average_Delta_0(
            np2vec(kR), 
            np2vec(E_kin_over_R), 
            H_delta_tau, 
            H_R, 
            expansion, 
            c_s,
            omega
        ); },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("H_delta_tau"), 
        py::arg("H_R"), 
        py::arg("expansion")=ExpansionHistory::RADIATION_DOMINATION,
        py::arg("c_s")=constants::cs_bag,
        py::arg("omega")=1.0/3.0,
        "Weighted average Delta0~ of Delta0 weighted with the kinetic spectrum"
        " Ekin.\n\n"
        "This is the weighted average of the low-frequency integral function"
        " Delta0 with the kinetic spectrum Ekin, parametrizing the dependence"
        " of the GW amplitude on the source duration delta_tau_fin, cf. Eq."
        " (63) of [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n"
        "    Momentum k times bubble separation R*\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k)/R*\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "H_R : float\n"
        "    Mean bubble separation in units of the Hubble size, H* R*\n"
        "expansion : ExpansionHistory, optional (default:"
        " ExpansionHistory.RADIATION_DOMINATION)\n"
        "    Assumed ExpansionHistory\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the stable phase\n"
        "omega : float, optional (default: 1/3)\n"
        "    Equation-of-state parameter omega of the Universe\n" 
        "\nReturns\n-------\n"
        "Delta0_tilde : float\n    Weighted average source lifetime dependence"
        " Delta0~(delta_tau_fin, R*)"
    );

    m.def(
        "normalized_GW_spectrum",
        [](
            py::array_t<double> Omega_GW, 
            double script_C, 
            double Omega_K,
            double Delta_tilde_0
        ) {
            return vec2np(normalized_GW_spectrum(
                np2vec(Omega_GW), 
                script_C, 
                Omega_K, 
                Delta_tilde_0
            ));
        },
        py::arg("Omega_GW"), 
        py::arg("script_C"), 
        py::arg("Omega_K"), 
        py::arg("Delta_tilde_0"), 
        "Normalized gravitational wave spectrum (k R*)^3 zeta_GW(k R*).\n\n"
        "This is the normalized gravitational wave spectrum in Eq. (93) of"
        " [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "Omega_GW : numpy.ndarray, shape (N,)\n"
        "    Gravitational wave power spectrum OmegaGW*(k) at production,"
        " cf. gravitational_wave_power_spectrum\n"
        "script_C : float\n"
        "    Amplitude parameter script_C of the anisotropic stress power"
        " spectrum, cf. anisotropic_stress_power_spectrum_amplitude\n"
        "Omega_K : float\n    Total kinetic energy density parameter Omega_K\n"
        "Delta_tilde_0 : float\n    Weighted average source lifetime dependence"
        " Delta0~, cf. weighted_average_Delta_0\n"
        "\nReturns\n-------\n"
        "K3_zeta_GW : numpy.ndarray, shape (N,)\n"
        "    Normalized gravitational wave spectrum K^3 zeta_GW"
    );

    m.def(
        "spectral_modification_function",
        [](py::array_t<double> zeta_GW, py::array_t<double> zeta_Pi) { 
            return vec2np(spectral_modification_function(
                np2vec(zeta_GW), np2vec(zeta_Pi)
            )); 
        },
        py::arg("zeta_GW"),  py::arg("zeta_Pi"),
        "Spectral modification function Delta~ = zeta_GW/zeta_Pi.\n\n"
        "This function parametrizes the modification of the spectral shape of"
        " the gravitational wave spectrum due to the finite source duration,"
        " Eq. (95) of [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "zeta_GW : numpy.ndarray, shape (N,)\n    Normalized gravitational wave"
        " spectrum zeta_GW (or K^3 zeta_GW) [cf. normalized_GW_spectrum]\n"
        "zeta_Pi : numpy.ndarray, shape (N,)\n"
        "    Normalized anisotropic stress spectrum zeta_Pi (or K^3 zeta_Pi)"
        " [cf. anisotropic_stress_power_spectral_shape]\n"
        "\nReturns\n-------\n"
        "Delta_tilde : numpy.ndarray, shape (N,)\n"
        "    Spectral modification function Delta~ = zeta_GW/zeta_Pi"
    );

    m.def(
        "upsilon_rad",
        &upsilon_rad,
        py::arg("H_delta_tau"),
        "Finite source lifetime suppression factor Upsilon in radiation"
        " domination.\n\n"
        "This is Eq. (85) of [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Source duration in units of the Hubble time,"
        " H* delta_tau_fin\n"
        "\nReturns\n-------\n"
        "Upsilon : float\n    Finite source lifetime suppression factor Upsilon"
    );

    m.def(
        "upsilon_flat",
        &upsilon_flat,
        py::arg("H_delta_tau"),
        "Finite source lifetime suppression factor Upsilon in a flat universe."
        "\n\nThis simply returns H* delta_tau_fin.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Source duration in units of the Hubble time,"
        " H* delta_tau_fin\n"
        "\nReturns\n-------\n"
        "Upsilon : float\n    Finite source lifetime suppression factor Upsilon"
    );

    m.def(
        "upsilon_const_eos",
        &upsilon_const_eos,
        py::arg("H_delta_tau"),
        py::arg("omega"),
        "Finite source lifetime suppression factor Upsilon in a universe with a"
        " constant equation of state.\n\n"
        "Cf. Eq. (3.12) of [Giombi+ (2025)].\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Source duration in units of the Hubble time,"
        " H* delta_tau_fin\n"
        "omega : float\n    Equation-of-state parameter omega of the Universe\n"
        "\nReturns\n-------\n"
        "Upsilon : float\n    Finite source lifetime suppression factor Upsilon"
    );

    m.def(
        "get_N_pR_default", 
        []() { return N_pR_default; }, 
        "Get the default number of points in the integration over p R*."
    );
    m.def(
        "set_N_pR_default", 
        [](size_t N) { N_pR_default = N; }, 
        py::arg("N"), 
        "Set the default number of points in the integration over p R*."
    );
    m.def(
        "get_N_z_default", 
        []() { return N_z_default; }, 
        "Get the default number of points in the integration over z."
    );
    m.def(
        "set_N_z_default", 
        [](size_t N) { N_z_default = N; }, 
        py::arg("N"), 
        "Set the default number of points in the integration over z."
    );
}