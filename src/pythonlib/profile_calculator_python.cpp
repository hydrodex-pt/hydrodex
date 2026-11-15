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
#include "hydrodex/profile_calculator.hpp"

namespace py = pybind11;
using namespace hydrodex;
using namespace hydrodex::profile_calculator;

/* Define macros to automate binding
 * Use recover_dynamic_default to control the default values with the
 * corresponding parameters of profile_calculator.
 */
// Macro for binding the alpha_min/max_* functions.
#define BIND_ALPHA_MIN_MAX(name, doc)                                          \
    m.def(#name,                                                               \
        py::vectorize([](                                                      \
            double xi_w,                                                       \
            double cp,                                                         \
            double cm,                                                         \
            py::object return_alpha_N,                                         \
            py::object step_size,                                              \
            py::object solve_mode                                              \
        ) {                                                                    \
            return name(                                                       \
                xi_w,                                                          \
                cp,                                                            \
                cm,                                                            \
                return_alpha_N.cast<bool>(),                                   \
                recover_dynamic_default(step_size, profile_step_size),         \
                recover_dynamic_default(solve_mode, profile_solve_mode)        \
            );                                                                 \
        }),                                                                    \
        py::arg("xi_w"),                                                       \
        py::arg("cp") = constants::cs_bag,                                     \
        py::arg("cm") = constants::cs_bag,                                     \
        py::arg("return_alpha_N") = true,                                      \
        py::arg("step_size") = py::none(),                                     \
        py::arg("solve_mode") = py::none(),                                    \
        doc                                                                    \
    )
// Macro for binding the compute_velocity_profile_* functions
#define BIND_COMPUTE_VELOCITY_PROFILE(name, arg1_name, arg2_name, doc)         \
    m.def(                                                                     \
        #name,                                                                 \
        [](                                                                    \
            double arg1_name,                                                  \
            double arg2_name,                                                  \
            double c_plus,                                                     \
            double c_minus,                                                    \
            py::object step_size,                                              \
            py::object solve_mode                                              \
        ) {                                                                    \
            std::vector<double> xi, v;                                         \
            size_t idx_wall = name(                                            \
                xi,                                                            \
                v,                                                             \
                arg1_name,                                                     \
                arg2_name,                                                     \
                c_plus,                                                        \
                c_minus,                                                       \
                recover_dynamic_default(step_size, profile_step_size),         \
                recover_dynamic_default(solve_mode, profile_solve_mode)        \
            );                                                                 \
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));            \
        },                                                                     \
        py::arg(#arg1_name),                                                   \
        py::arg(#arg2_name),                                                   \
        py::arg("c_plus") = constants::cs_bag,                                 \
        py::arg("c_minus") = constants::cs_bag,                                \
        py::arg("step_size") = py::none(),                                     \
        py::arg("solve_mode") = py::none(),                                    \
        doc                                                                    \
    )
// end of the macros

void bind_profile_calculator(pybind11::module_& m) {

    // Hydrodynamic solution types
    py::native_enum<HydroSolutionType>(
        m, 
        "HydroSolutionType", 
        "enum.Enum",
        "Types of the hydrodynamic solutions."
    )
        .value(
            "Deflagration", 
            HydroSolutionType::Deflagration,
            "weak (subsonic) direct deflagration"
        )
        .value(
            "Detonation",
            HydroSolutionType::Detonation,
            "weak (supersonic) direct detonation"
        )
        .value(
            "Hybrid",
            HydroSolutionType::Hybrid,
            "direct hybrid solution (supersonic deflagration, '(-)' branch)"
        )
        .value(
            "HybridPlus",
            HydroSolutionType::HybridPlus,
            "direct hybrid solution (supersonic deflagration, '(+)' branch)"
        )
        .value(
            "InverseDeflagration",
            HydroSolutionType::InverseDeflagration,
            "weak (supersonic) inverse deflagration"
        )
        .value(
            "InverseDetonation",
            HydroSolutionType::InverseDetonation,
            "weak (subsonic) inverse detonation ('(-)' branch)"
        )
        .value(
            "InverseDetonationPlus",
            HydroSolutionType::InverseDetonationPlus,
            "weak (subsonic) inverse detonation ('(+)' branch)"
        )
        .value(
            "InverseHybrid",
            HydroSolutionType::InverseHybrid,
            "inverse hybrid solution (subsonic inverse deflagration)"
        )
        .value(
            "MixedHybridMinus",
            HydroSolutionType::MixedHybridMinus,
            "mixed direct-inverse hybrid solution (v- = c-)"
        )
        .value(
            "MixedHybridPlus",
            HydroSolutionType::MixedHybridPlus,
            "mixed direct-inverse hybrid solution (v+ = c+)"
        )
    .finalize();
    auto hydrosol_type = m.attr("HydroSolutionType");
    hydrosol_type.attr("to_string") = py::cpp_function(
        [](py::object self) {return to_string(self.cast<HydroSolutionType>());},
        py::is_method(hydrosol_type)
    );

    // Efficiency-factor normalization conventions
    py::native_enum<EfficiencyNormalization>(
        m,
        "EfficiencyNormalization",
        "enum.Enum",
        "Normalization conventions for the hydrodynamic efficiency factor."
    )
        .value(
            "BubbleVolume",
            EfficiencyNormalization::BubbleVolume,
            "Normalize to the bubble volume. This corresponds to the kinetic"
            " energy density rho_kin."
        )
        .value(
            "EnergyDensity",
            EfficiencyNormalization::TotalEnergyDensity,
            "Normalize to the total energy density. This corresponds to the"
            " energy fraction K = rho_kin/rho_tot. (Deprecated name: use"
            " TotalEnergyDensity, which matches the C++ enumerator.)"
        )
        .value(
            "TotalEnergyDensity",
            EfficiencyNormalization::TotalEnergyDensity,
            "Normalize to the total energy density. This corresponds to the"
            " energy fraction K = rho_kin/rho_tot (same as EnergyDensity)."
        )
        .value(
            "VacuumEnergy",
            EfficiencyNormalization::VacuumEnergy,
            "Normalize to the vacuum energy density. This corresponds to"
            " kappa = rho_kin/rho_vac with rho_vac = 4/3 pi xi_w^3 epsilon. For"
            " negative epsilon < 0, this option normalizes to the initial"
            " thermal energy rho_ini = 4/3 pi xi_w^3 wN/(1+c+^2) instead."
        )
        .value(
            "AlphaPseudoTrace",
            EfficiencyNormalization::AlphaPseudoTrace,
            "Normalize to the pseudotrace strength parameter alphaN at"
            " nucleation. This corresponds to kappa_bar = 4 rho_kin / (3 alphaN"
            " wN). For negative alphaN < 0, this option normalizes to the"
            " initial thermal energy rho_ini = 4/3 pi xi_w^3 wN/(1+c+^2)"
            " instead."
        )
    .finalize();

    // Profile calculation modes 
    py::native_enum<ProfileSolveMode>(
        m, 
        "ProfileSolveMode", 
        "enum.Enum", 
        "Selection which differential equation to solve when calculating the"
        " velocity profile."
    )
        .value("D_XI_D_V",  
            ProfileSolveMode::D_XI_D_V,  
            "always solve d(xi)/d(v)"
        )
        .value(
            "D_V_D_XI",  
            ProfileSolveMode::D_V_D_XI,  
            "always solve d(v)/d(xi)"
        )
        .value(
            "D_V_D_XI_FOR_SHOCKS",  
            ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,  
            "solve d(v)/d(xi) when the solution is expected to end in a shock"
            " and d(xi)/d(v) otherwise"
        )
    .finalize();

    // Root branch selector in velocity matching conditions
    py::native_enum<RootBranch>(
        m, 
        "RootBranch", 
        "enum.Enum", 
        "Branch selector for the sign in front of the square root in the"
        " matching for v+ and v-."
    )
        .value(
            "Plus", 
            RootBranch::Plus, 
            "branch with positive sign in front of the square root"
        )
        .value(
            "Minus", 
            RootBranch::Minus, 
            "branch with negative sign in front of the square root"
        )
    .finalize();

    m.def(
        "gamma",
        py::vectorize(&profile_calculator::gamma),
        py::arg("v"),
        "Compute the Lorentz boost factor gamma=1/sqrt(1-v^2).\n"
        "\nParameters\n----------\n"
        "v : float or array_like\n    Velocity v\n"
        "\nReturns\n-------\n"
        "gamma : float or ndarray\n    Boost factor gamma"
    );

    m.def(
        "mu",
        py::vectorize(&profile_calculator::mu),
        py::arg("xi"), py::arg("v"),
        "Lorentz boosted velocity.\n\n"
        "This boosts the velocity v to the frame moving with velocity xi,"
        " Eq. (2.7) of [Barni+ (2027)].\n"
        "    mu = (xi - v)/(1 - xi v)\n"
        "\nParameters\n----------\n"
        "xi : float or array_like\n    Relative velocity xi of the new frame"
        " to which the velocity is boosted\n"
        "v : float or array_like\n    Velocity v in the original frame\n"
        "\nReturns\n-------\n"
        "mu : float or ndarray\n    Boosted velocity mu(xi,v)"
    );

    m.def(
        "check_velocity",
        &check_velocity,
        py::arg("v"),
        "Check if the velocity is between 0 and 1.\n"
        "\nParameters\n----------\n"
        "v : float\n    Velocity to check\n"
        "\nReturns\n-------\n"
        "ok : bool\n"
        "    True if the velocity lies between 0 and 1, False otherwise"
    );

    m.def(
        "get_nu_from_cs",
        &get_nu_from_cs,
        py::arg("cs"),
        "Get the temperature exponent nu from the speed of sound cs in the"
        " constant-sound-speed model.\n"
        "\nParameters\n----------\n"
        "cs : float\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "nu : float\n    Temperature exponent nu"
    );

    m.def(
        "get_cs_from_nu",
        &get_cs_from_nu,
        py::arg("nu"),
        "Get the speed of sound cs from the temperature exponent nu in the"
        " constant-sound-speed model.\n"
        "\nParameters\n----------\n"
        "nu : float\n    Temperature exponent nu\n"
        "\nReturns\n-------\n"
        "cs : float\n    Speed of sound cs"
    );

    m.def(
        "shock_condition",
        &shock_condition,
        py::arg("xi"), py::arg("v"), py::arg("cs")=constants::cs_bag,
        "Calculate the distance to the shock front.\n\n"
        "This function returns the shock condition mu(xi_sh, v_sh) xi_sh - cs^2"
        ". It returns zero if (xi,v) corresponds to a shock front, a negative"
        " value if the point is ahead of a shock front and a positive value if"
        " it is behind the shock front.\n"
        "\nParameters\n----------\n"
        "xi : float\n    Self-similar coordinate xi = r/t\n"
        "v : float\n    Fluid velocity v\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "delta : float\n    Distance to the shock front"
    );

    m.def(
        "jouguet_condition",
        &jouguet_condition,
        py::arg("xi"), py::arg("v"), py::arg("cs")=constants::cs_bag,
        "Calculate the distance to the Jouguet front.\n\n"
        "This function returns the Jouguet condition mu(xi_J, v_J)^2 - cs^2."
        " It returns zero if (xi,v) corresponds to a Jouguet front, a negative"
        " value if the point is ahead of a Jouguet front and a positive value"
        " if it is behind the Jouguet front.\n"
        "\nParameters\n----------\n"
        "xi : float\n    Self-similar coordinate xi = r/t\n"
        "v : float\n    Fluid velocity v\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "delta : float\n    Distance to the Jouguet front"
    );

    m.def(
        "get_v_plus_from_v_minus",
        [](double vm, double alpha, double cm) { 
            auto [v1, v2] = get_v_plus_from_v_minus(vm, alpha, cm);
            return py::make_tuple(v1,v2); 
        },
        py::arg("vm"), py::arg("alpha"), py::arg("cm")=constants::cs_bag,
        "Calculate v+ from v- and alpha+.\n\n"
        "This function calculates the fluid velocity v+ ahead of the bubble (or"
        " shock) front from the velocity v- behind the front and the"
        " pseudotrace transition strength parameter alpha+ ahead of the"
        " front.\n\n"
        "This function returns both branches of the solution for v+.\n"
        "\nParameters\n----------\n"
        "vm : float\n    Fluid velocity v- behind the front\n"
        "alpha : float\n    "
        "Pseudotrace transition strength parameter alpha+ ahead of the front\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the front\n"
        "\nReturns\n-------\n"
        "vp1 : float\n    The solution on the '(-)' branch for the fluid"
        " velocity v+ ahead of the front\n"
        "vp2 : float\n    The solution on the '(+)' branch for the fluid"
        " velocity v+ ahead of the front"
    );

    m.def(
        "get_v_plus_from_v_minus",
        py::vectorize(
            py::overload_cast<double, double, RootBranch, double>(
                &get_v_plus_from_v_minus
            )
        ),
        py::arg("vm"),
        py::arg("alpha"),
        py::arg("branch"),
        py::arg("cm")=constants::cs_bag,
        "Calculate v+ from v- and alpha+.\n\n"
        "This function calculates the fluid velocity v+ ahead of the bubble (or"
        " shock) front from the velocity v- behind the front and the"
        " pseudotrace transition strength parameter alpha+ ahead of the"
        " front.\n\n"
        "This function returns the '(-)' or '(+)' branch of the solution for"
        " v+ depending on the value of the `branch` argument.\n"
        "\nParameters\n----------\n"
        "vm : float or array_like\n    Fluid velocity v- behind the front\n"
        "alpha : float or array_like\n    "
        "Pseudotrace transition strength parameter alpha+ ahead of the front\n"
        "branch : RootBranch\n"
        "    RootBranch selector to determine which branch to return\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the front\n"
        "\nReturns\n-------\n"
        "vp : float or ndarray\n    Solution for the fluid velocity v+ ahead"
        " of the front on the selected branch"
    );

    m.def(
        "get_v_minus_from_v_plus",
        [](double vp, double alpha, double cm) { 
            auto [v1, v2] = get_v_minus_from_v_plus(vp, alpha, cm);
            return py::make_tuple(v1,v2); 
        },
        py::arg("vp"), py::arg("alpha"), py::arg("cm")=constants::cs_bag,
        "Calculate v- from v+ and alpha+.\n\n"
        "This function calculates the fluid velocity v- behind the bubble (or"
        " shock) front from the velocity v+ and the pseudotrace transition"
        " strength alpha+ ahead of the front.\n\n"
        "This function returns both branches of the solution for v-.\n"
        "\nParameters\n----------\n"
        "vp : float\n    Fluid velocity v+ ahead of the front\n"
        "alpha : float\n    "
        "Pseudotrace transition strength parameter alpha+ ahead of the front\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the front\n"
        "\nReturns\n-------\n"
        "vm1 : float\n    The solution on the '(-)' branch for the fluid"
        " velocity v- behind the front\n"
        "vm2 : float\n    The solution on the '(+)' branch for the fluid"
        " velocity v- behind the front"
    );

    m.def(
        "get_v_minus_from_v_plus",
        py::vectorize(
            py::overload_cast<double, double, RootBranch, double>(
                &get_v_minus_from_v_plus
            )
        ),
        py::arg("vp"),
        py::arg("alpha"),
        py::arg("branch"),
        py::arg("cm")=constants::cs_bag,
        "Calculate v- from v+ and alpha+.\n\n"
        "This function calculates the fluid velocity v- behind the bubble (or"
        " shock) front from the velocity v+ and the pseudotrace transition"
        " strength alpha+ ahead of the front.\n\n"
        "This function returns the '(-)' or '(+)' branch of the solution for"
        " v- depending on the value of the `branch` argument.\n"
        "\nParameters\n----------\n"
        "vp : float or array_like\n    Fluid velocity v+ ahead of the front\n"
        "alpha : float or array_like\n    "
        "Pseudotrace transition strength parameter alpha+ ahead of the front\n"
        "branch : RootBranch\n"
        "    RootBranch selector to determine which branch to return\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the front\n"
        "\nReturns\n-------\n"
        "vm : float or ndarray\n    Solution for the fluid velocity v- behind"
        " the front on the selected branch"
    );

    m.def(
        "compute_alpha_plus",
        py::vectorize(&compute_alpha_plus),
        py::arg("vp"), py::arg("vm"), py::arg("cm")=constants::cs_bag,
        "Compute alpha+ from v+ and v-.\n\n"
        "This function calculates the pseudotrace transition strength ahead of"
        " the wall, alpha+, from the fluid velocities v+ and v- in the wall"
        " rest frame.\n\n"
        "For c- = cs_bag = 1/sqrt(3), this coincides with the bag-model"
        " transition strength alpha_bag+.\n"
        "\nParameters\n----------\n"
        "vp : float or array_like\n    Fluid velocity v+ ahead of the wall\n"
        "vm : float or array_like\n    Fluid velocity v- behind the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "alpha : float or ndarray\n    Pseudotrace phase transition strength"
        " alpha+ ahead of the bubble wall"
    );

    m.def(
        "get_v_plus_extremum_from_v_minus",
        py::vectorize(&get_v_plus_extremum_from_v_minus),
        py::arg("vm"), py::arg("cm")=constants::cs_bag,
        "Get the extremum position v+(extr).\n\n"
        "The '(-)' and '(+)' branch of v-(v+, alpha+) exhibit a maximum or"
        " minimum, respectively, for strength parameters "
        "- (1-c-^2)/(6 c-^2) < alpha+ < 0.\n\n"
        "This function returns the position v+(extr) of the extremum, cf."
        " Eq. (3.9a) of [Barni+ (2027)], with the value of alpha+ such that"
        " the maximal or minimal value is given by the input parameter vm.\n"
        "\nParameters\n----------\n"
        "vm : float or array_like\n"
        "    Value v-(extr) of v- at the extremum\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "vp : float or ndarray\n    Position v+(extr) of the extremum"
    );

    m.def(
        "get_v_plus_extremum_from_alpha",
        py::vectorize(&get_v_plus_extremum_from_alpha),
        py::arg("alpha"), py::arg("cm")=constants::cs_bag,
        "Get the extremum position v+(extr).\n\n"
        "The '(-)' and '(+)' branch of v-(v+, alpha+) exhibit a maximum or"
        " minimum, respectively, for strength parameters "
        "- (1-c-^2)/(6 c-^2) < alpha+ < 0.\n\n"
        "This function returns the position v+(extr) of the extremum, cf."
        " Eq. (3.9a) of [Barni+ (2027)], for a given value of alpha+.\n"
        "\nParameters\n----------\n"
        "alpha : float or array_like\n    "
        "Pseudotrace transition strength parameter alpha+ ahead of the front\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "vp : float or ndarray\n    Position v+(extr) of the extremum"
    );

    m.def(
        "get_v_minus_extremum_from_v_plus",
        py::vectorize(&get_v_minus_extremum_from_v_plus),
        py::arg("vp"), py::arg("branch"), py::arg("cm")=constants::cs_bag,
        "Get the extremum value v-(extr).\n\n"
        "The '(-)' and '(+)' branch of v-(v+, alpha+) exhibit a maximum or"
        " minimum, respectively, for strength parameters "
        "- (1-c-^2)/(6 c-^2) < alpha+ < 0.\n\n"
        "This function returns the value v-(extr) at the extremum on the"
        " branch determined by `branch`, cf. Eq. (3.9b) of [Barni+ (2027)],"
        " with alpha+ such that the extremum occurs at the position v+(extr)"
        " given by the input parameter vp.\n"
        "\nParameters\n----------\n"
        "vp : float or array_like\n    Position v+(extr) of the extremum\n"
        "branch : RootBranch\n"
        "    RootBranch selector to determine which branch to return\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "vm : float or ndarray\n    Value v-(extr) at the extremum"
    );

    m.def(
        "get_v_minus_extremum_from_alpha",
        py::vectorize(&get_v_minus_extremum_from_alpha),
        py::arg("alpha"), py::arg("branch"), py::arg("cm")=constants::cs_bag,
        "Get the extremum value v-(extr).\n\n"
        "The '(-)' and '(+)' branch of v-(v+, alpha+) exhibit a maximum or"
        " minimum, respectively, for strength parameters "
        "- (1-c-^2)/(6 c-^2) < alpha+ < 0.\n\n"
        "This function returns the value v-(extr) at the extremum on the"
        " branch determined by `branch`, cf. Eq. (3.9b) of [Barni+ (2027)],"
        " for a given value of alpha+.\n"
        "\nParameters\n----------\n"
        "alpha : float or array_like\n    "
        "Pseudotrace transition strength parameter alpha+ ahead of the front\n"
        "branch : RootBranch\n"
        "    RootBranch selector to determine which branch to return\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "vm : float or ndarray\n    Value v-(extr) at the extremum"
    );

    m.def(
        "direct_jouguet_velocity",
        py::vectorize(&direct_jouguet_velocity),
        py::arg("alpha"),
        py::arg("cm")=constants::cs_bag,
        "Compute the Jouguet velocity for a direct phase transition.\n"
        "\nParameters\n----------\n"
        "alpha : float or array_like\n"
        "    Pseudotrace strength parameter alpha+\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n" 
        "vJ : float or ndarray\n    Jouguet velocity vJ(direct)"
    );

    m.def(
        "inverse_jouguet_velocity",
        py::vectorize(&inverse_jouguet_velocity),
        py::arg("alpha"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Compute the Jouguet velocity for an inverse phase transition.\n"
        "\nParameters\n----------\n"
        "alpha : float or array_like\n"
        "    Pseudotrace strength parameter alpha+\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n" 
        "vJinv : float or ndarray\n    Inverse Jouguet velocity vJ(inverse)"
    );

    m.def(
        "xi_min",
        py::vectorize(&xi_min),
        py::arg("sol_type"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Minimal wall velocity xi_w for a given solution type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n" 
        "xi_min : float or ndarray\n"
        "    The minimal wall velocity for the solution type"
    );

    m.def(
        "xi_max",
        py::vectorize(&xi_max),
        py::arg("sol_type"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Maximal wall velocity xi_w for a given solution type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n" 
        "xi_max : float or ndarray\n"
        "    The maximal wall velocity for the solution type"
    );

    m.def(
        "check_wall_velocity",
        &check_wall_velocity,
        py::arg("sol_type"),
        py::arg("xi_w"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Check if `xi_w` is a valid wall velocity for solutions of type"
        " `sol_type`.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nRaises\n------\n"
        "ValueError\n    If the wall velocity is not valid"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_detonation,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for detonation solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct detonations.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " direct detonations"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_detonation,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for detonation solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct detonations.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " direct detonations"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_hybrid_plus,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for direct hybrid ('(+)' branch) solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct hybrids on the '(+)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " direct hybrids ('(+)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_hybrid_plus,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for direct hybrid ('(+)' branch) solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct hybrids on the '(+)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " direct hybrids ('(+)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_hybrid,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for direct hybrid ('(-)' branch) solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct hybrids on the '(-)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " direct hybrids ('(-)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_hybrid,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for direct hybrid ('(-)' branch) solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct hybrids on the '(-)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " direct hybrids ('(-)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_deflagration,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for deflagration solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct deflagrations.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " direct deflagrations"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_deflagration,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for deflagration solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for direct deflagrations.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " direct deflagrations"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_inverse_detonation,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse detonation ('(-)' branch) solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse detonations on the '(-)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " inverse detonations ('(-)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_inverse_detonation,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse detonation ('(-)' branch) solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse detonations on the '(-)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " inverse detonations ('(-)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_inverse_detonation_plus,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse detonation ('(+)' branch) solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse detonations on the '(+)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " inverse detonations ('(+)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_inverse_detonation_plus,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse detonation ('(+)' branch) solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse detonations on the '(+)' branch.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " inverse detonations ('(+)' branch)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_inverse_hybrid,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse hybrid solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse hybrids.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " inverse hybrids"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_inverse_hybrid,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse hybrid solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse hybrids.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " inverse hybrids"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_inverse_deflagration,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse deflagration solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse deflagrations.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " inverse deflagrations"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_inverse_deflagration,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for inverse deflagration solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for inverse deflagrations.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " inverse deflagrations"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_mixed_hybrid_minus,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for mixed hybrid (v- = c-) solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for mixed hybrids with v- = c-.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " mixed hybrids (v- = c-)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_mixed_hybrid_minus,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for mixed hybrid (v- = c-) solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for mixed hybrids with v- = c-.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " mixed hybrids (v- = c-)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_max_mixed_hybrid_plus,
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for mixed hybrid (v+ = c+) solutions.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for mixed hybrids with v+ = c+.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " mixed hybrids (v+ = c+)"
    );

    BIND_ALPHA_MIN_MAX(
        alpha_min_mixed_hybrid_plus,
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for mixed hybrid (v+ = c+) solutions.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for mixed hybrids with v+ = c+.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " mixed hybrids (v+ = c+)"
    );

    m.def(
        "alpha_max",
        py::vectorize([](
            HydroSolutionType sol_type, 
            double xi_w,
            double cp, 
            double cm, 
            py::object return_alpha_N, 
            py::object step_size, 
            py::object solve_mode
        ) {
            return alpha_max(
                sol_type,
                xi_w,
                cp,
                cm,
                return_alpha_N.cast<bool>(),
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        }),
        py::arg("sol_type"), 
        py::arg("xi_w"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+ for the given solution type.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for the solution type specified in sol_type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength for"
        " the specified solution type"
    );

    m.def(
        "alpha_min",
        py::vectorize([](
            HydroSolutionType sol_type, 
            double xi_w,
            double cp, 
            double cm, 
            py::object return_alpha_N, 
            py::object step_size, 
            py::object solve_mode
        ) {
            return alpha_min(
                sol_type,
                xi_w,
                cp,
                cm,
                return_alpha_N.cast<bool>(),
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        }),
        py::arg("sol_type"), 
        py::arg("xi_w"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+ for the given solution type.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) for the solution type specified in sol_type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength for"
        " the specified solution type"
    );

    m.def(
        "alpha_max",
        py::vectorize([](
            double xi_w,
            double cp, 
            double cm, 
            py::object return_alpha_N, 
            py::object step_size, 
            py::object solve_mode
        ) {
            return alpha_max(
                xi_w,
                cp,
                cm,
                return_alpha_N.cast<bool>(),
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        }),
        py::arg("xi_w"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute the maximal value of the pseudotrace transition strength"
        " alphaN or alpha+.\n\n"
        "This function returns the maximal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) over all solution types for the given wall velocity.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float or ndarray\n"
        "    The maximal value of the transition strength"
    );

    m.def(
        "alpha_min",
        py::vectorize([](
            double xi_w,
            double cp, 
            double cm, 
            py::object return_alpha_N, 
            py::object step_size, 
            py::object solve_mode
        ) {
            return alpha_min(
                xi_w,
                cp,
                cm,
                return_alpha_N.cast<bool>(),
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        }),
        py::arg("xi_w"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute the minimal value of the pseudotrace transition strength"
        " alphaN or alpha+.\n\n"
        "This function returns the minimal value of the pseudotrace phase"
        " transition strength alphaN far ahead of the wall (if return_alpha_N"
        " is True) or alpha+ directly in front of the wall (if return_alpha_N"
        " is False) over all solution types for the given wall velocity.\n"
        "\nParameters\n----------\n"
        "xi_w : float or array_like\n    Wall velocity xi_w\n"
        "cp : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return min(alphaN) (True) or min(alpha+) (False)\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float or ndarray\n"
        "    The minimal value of the transition strength"
    );

    m.def(
        "inverse_detonation_plus_fold_critical_cm",
        py::vectorize(&inverse_detonation_plus_fold_critical_cm),
        py::arg("cp"),
        "Critical value of c- for onset of non-injectivity on the '(+)' branch"
        " of inverse detonations.\n\n"
        "Returns the fitted critical sound speed\n\n"
        "    cm_crit(cp) = cp [1 - (1 - cp^2) P(cp^2)],\n\n"
        "where P(x) is a degree-20 polynomial fit.\n"
        "\nParameters\n----------\n"
        "cp : float or array_like\n    Speed of sound c+ ahead of the wall\n"
        "\nReturns\n-------\n"
        "cm_crit : float or ndarray\n"
        "    Critical value cm_crit(cp) at which non-injectivity sets in"
    );

    m.def(
        "check_noninjectivity",
        py::vectorize(&check_noninjectivity),
        py::arg("cp"),
        py::arg("cm"),
        "Check whether (c+, c-) lies in the non-injective sound-speed region.\n"
        "\nReturns True if c+ > c- and c- < cm_crit(c+), where cm_crit(c+) is"
        " the fitted boundary from inverse_detonation_plus_fold_critical_cm.\n"
        "\nParameters\n----------\n"
        "cp : float or array_like\n    Speed of sound c+ ahead of the wall\n"
        "cm : float or array_like\n    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "noninjective : bool or ndarray\n    True if the sound speeds lie in"
        " the region where more than one solution can exist for some"
        " (xi_w, alphaN) values"
    );

    m.def(
        "classify_wave_fluid_frame",
        [](
            double xi_wall,
            double alpha,
            double cp, 
            double cm, 
            py::object step_size, 
            py::object solve_mode
        ) {
            return classify_wave_fluid_frame(
                xi_wall,
                alpha,
                cp,
                cm,
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        },
        py::arg("xi_wall"), 
        py::arg("alpha"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Classify the type of the solution based on xi_wall and the pseudotrace"
        " strength parameter alphaN in the fluid frame.\n\n"
        "Unless the noninjectivity check mode (see set_noninjectivity_check) is"
        " set to CheckMode.SILENT, the sound-speed pair is checked with"
        " profile_calculator.check_noninjectivity(cp, cm) before branch"
        " selection and a warning or error is thrown depending on the"
        " noninjectivity check mode.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Bubble wall velocity xi_w\n"
        "alpha : float\n"
        "    Pseudotrace phase transition strength alphaN (far ahead of the"
        " wall)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "hydrotype : HydroSolutionType\n    Type of the hydrodynamic solution"
    );

    m.def(
        "classify_wave_wall_frame",
        &classify_wave_wall_frame,
        py::arg("vp"), 
        py::arg("vm"), 
        py::arg("xi_w"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag, 
        "Classify the type of the solution based on v+ and v- in the wall rest"
        " frame and the wall velocity xi_w.\n"
        "\nParameters\n----------\n"
        "vp : float\n    Fluid velocity v+ ahead of the wall\n"
        "vm : float\n    Fluid velocity v- behind the wall\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "hydrotype : HydroSolutionType\n    Type of the hydrodynamic solution"
    );

    m.def(
        "compute_fluid_velocities",
        [](
            double xi_wall, 
            double alpha, 
            HydroSolutionType sol_type, 
            double cp, 
            double cm
        ) { 
            auto [vp, vm] = compute_fluid_velocities(
                xi_wall, alpha, sol_type, cp, cm
            );
            return py::make_tuple(vp,vm); 
        },
        py::arg("xi_wall"), 
        py::arg("alpha"),
        py::arg("sol_type"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag, 
        "Compute the fluid velocities v+ and v- in the wall rest frame.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Bubble wall velocity xi_w\n"
        "alpha : float\n"
        "    Pseudotrace phase transition strength alpha+ ahead of the wall\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "vp : float\n"
        "    The fluid velocity v+ ahead of the wall in the wall rest frame\n"
        "vm : float\n"
        "    The fluid velocity v- behind the wall in the wall rest frame"
    );

    m.def(
        "dxi_dv",
        &dxi_dv,
        py::arg("v"), py::arg("xi"), py::arg("cs")=constants::cs_bag,
        "Derivative dxi/dv.\n\n"
        "We solve this differential equation to obtain xi(v) and then invert to"
        " obtain v(xi). This is Eq. (2.10) of [Barni+ (2027)].\n"
        "\nParameters\n----------\n"
        "v : float\n    Fluid velocity v\n"
        "xi : float\n    Self-similar coordinate xi=r/t\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "dxi_dv : float\n    Derivative d(xi)/d(v)"
    );

    m.def(
        "dv_dxi",
        &dv_dxi,
        py::arg("xi"), py::arg("v"), py::arg("cs")=constants::cs_bag,
        "Derivative dv/dxi.\n\n"
        "We solve this differential equation to directly obtain v(xi)."
        " This is Eq. (2.10) of [Barni+ (2027)].\n"
        "\nParameters\n----------\n"
        "xi : float\n    Self-similar coordinate xi=r/t\n"
        "v : float\n    Fluid velocity v\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "dv_dxi : float\n    Derivative d(v)/d(xi)"
    );

    m.def(
        "integrate_xi_of_v",
        [](
            py::array_t<double> xi_arr, 
            py::array_t<double> v_arr, 
            double v_stop, 
            double dv, 
            bool stop_at_shock,
            double cs
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = integrate_xi_of_v(
                xi_vec, v_vec, v_stop, dv, stop_at_shock, cs
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"), 
        py::arg("v_arr"), 
        py::arg("v_stop"), 
        py::arg("dv"), 
        py::arg("stop_at_shock")=true,
        py::arg("cs")=constants::cs_bag,
        "Solve the differential equation d(xi)/d(v) for the velocity profile."
        "\n\n"
        "This function solves the differential equation using a fourth-order"
        " Runge-Kutta method. The first two arguments are arrays to which the"
        " solutions are appended. These arrays must already contain the initial"
        " condition.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "v_stop : float\n    Final fluid velocity vf\n"
        "dv : float\n    Fluid velocity step size\n"
        "stop_at_shock : bool, optional (default: True)\n"
        "    Whether to stop the evolution when a shock front is encountered\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    A status code indicating the reason for termination:\n"
        "        0  - reached v_stop\n"
        "        1  - encountered shock\n"
        "        -1 - entered unphysical region for xi\n"
        "        -2 - reached superluminal fluid velocity v\n"
        "xi_arr : numpy.ndarray, shape (N',)\n  "
        "  Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended"
    );

    m.def(
        "integrate_v_of_xi",
        [](
            py::array_t<double> xi_arr, 
            py::array_t<double> v_arr, 
            double v_stop, 
            double dxi, 
            bool stop_at_shock,
            double cs
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = integrate_v_of_xi(
                xi_vec, v_vec, v_stop, dxi, stop_at_shock, cs
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"), 
        py::arg("v_arr"), 
        py::arg("v_stop"), 
        py::arg("dxi"), 
        py::arg("stop_at_shock")=true,
        py::arg("cs")=constants::cs_bag,
        "Solve the differential equation d(v)/d(xi) for the velocity profile."
        "\n\n"
        "This function solves the differential equation using a fourth-order"
        " Runge-Kutta method. The first two arguments are arrays to which the"
        " solutions are appended. These arrays must already contain the initial"
        " condition.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "v_stop : float\n    Final fluid velocity vf\n"
        "dxi : float\n    Self-similar coordinate step size\n"
        "stop_at_shock : bool, optional (default: True)\n"
        "    Whether to stop the evolution when a shock front is encountered\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    A status code indicating the reason for termination:\n"
        "        0  - reached v_stop\n"
        "        1  - encountered shock\n"
        "        -1 - entered unphysical region for xi\n"
        "        -2 - reached superluminal fluid velocity v\n"
        "        -3 - crossed Jouguet front\n"
        "xi_arr : numpy.ndarray, shape (N',)\n  "
        "  Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended"
    );

    m.def(
        "integrate_velocity_profile",
        [](
            py::array_t<double> xi_arr, 
            py::array_t<double> v_arr, 
            double v_stop, 
            double step_size, 
            bool stop_at_shock,
            double cs,
            py::object solve_mode
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = integrate_velocity_profile(
                xi_vec, v_vec, v_stop, step_size, stop_at_shock, cs, 
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"), 
        py::arg("v_arr"), 
        py::arg("v_stop"), 
        py::arg("step_size"), 
        py::arg("stop_at_shock")=true,
        py::arg("cs")=constants::cs_bag,
        py::arg("solve_mode")=py::none(),
        "Solve the differential equation for the velocity profile.\n\n"
        "This function solves the differential equation using a fourth-order"
        " Runge-Kutta method. The first two arguments are arrays to which the"
        " solutions are appended. These arrays must already contain the initial"
        " condition.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "v_stop : float\n    Final fluid velocity vf\n"
        "step_size : float\n    Fluid velocity or self-similar coordinate step"
        " size (depending on solve_mode)\n"
        "stop_at_shock : bool, optional (default: True)\n"
        "    Whether to stop the evolution when a shock front is encountered\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    A status code indicating the reason for termination:\n"
        "        0  - reached v_stop\n"
        "        1  - encountered shock\n"
        "        -1 - entered unphysical region for xi\n"
        "        -2 - reached superluminal fluid velocity v\n"
        "        -3 - crossed Jouguet front\n"
        "xi_arr : numpy.ndarray, shape (N',)\n  "
        "  Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended"
    );

    m.def(
        "compute_rarefaction_wave",
        [](
            py::array_t<double> xi_arr,
            py::array_t<double> v_arr,
            double c_s,
            py::object step_size,
            py::object solve_mode
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = compute_rarefaction_wave(
                xi_vec, v_vec, c_s,
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"),
        py::arg("v_arr"),
        py::arg("c_s")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute a rarefaction wave attached to the bubble wall.\n\n"
        "This function integrates the rarefaction-wave part of a velocity"
        " profile from the bubble wall, given by the last elements of xi_arr and"
        " v_arr, towards the speed of sound xi = c_s.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi (must contain the"
        " wall position xi_w as last element)\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v (must contain the fluid"
        " velocity v_w at the wall as last element)\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the region of the rarefaction wave\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    Status code returned by integrate_velocity_profile (0 if the"
        " integration reached the final fluid velocity)\n"
        "xi_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended\n"
        "\nRaises\n------\n"
        "RuntimeError\n"
        "    If the initial conditions are not valid for a rarefaction wave"
    );

    m.def(
        "compute_compression_wave",
        [](
            py::array_t<double> xi_arr,
            py::array_t<double> v_arr,
            double c_s,
            py::object step_size,
            py::object solve_mode
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = compute_compression_wave(
                xi_vec, v_vec, c_s,
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"),
        py::arg("v_arr"),
        py::arg("c_s")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute a compression wave attached to the bubble wall.\n\n"
        "This function integrates the compression-wave part of a velocity"
        " profile from the bubble wall, given by the last elements of xi_arr and"
        " v_arr, until a shock front is encountered.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi (must contain the"
        " wall position xi_w as last element)\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v (must contain the fluid"
        " velocity v_w at the wall as last element)\n"
        "c_s : float, optional (default: cs_bag)\n"
        "    Speed of sound c_s in the region of the compression wave\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    Status code returned by integrate_velocity_profile (1 if a shock"
        " front was encountered)\n"
        "xi_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended\n"
        "\nRaises\n------\n"
        "RuntimeError\n"
        "    If the initial conditions are not valid for a compression wave"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_direct_detonation,
        v_plus,
        v_minus,
        "Compute the velocity profile for a direct detonation.\n\n"
        "This function solves the fluid velocity profile for a direct"
        " detonation and returns the results of the velocity v as a function of"
        " the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n" 
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_direct_deflagration,
        v_plus, 
        v_minus, 
        "Compute the velocity profile for a direct deflagration.\n\n"
        "This function solves the fluid velocity profile for a direct"
        " deflagration and returns the results of the velocity v as a function"
        " of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_direct_hybrid,
        v_plus,
        xi_wall,
        "Compute the velocity profile for a direct hybrid ('(+/-)' branch).\n\n"
        "This function solves the fluid velocity profile for direct hybrids"
        " (on either branch) and returns the results of the velocity v as a"
        " function of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_inverse_detonation,
        v_plus,
        v_minus,
        "Compute the velocity profile for an inverse detonation ('(+/-)'"
        " branch).\n\n"
        "This function solves the fluid velocity profile for an inverse"
        " detonation on either branch and returns the results of the velocity"
        " v as a function of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_inverse_deflagration,
        v_plus,
        v_minus,
        "Compute the velocity profile for an inverse deflagration.\n\n"
        "This function solves the fluid velocity profile for an inverse"
        " deflagration and returns the results of the velocity v as a function"
        " of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_inverse_hybrid,
        xi_wall,
        v_minus,
        "Compute the velocity profile for an inverse hybrid.\n\n"
        "This function solves the fluid velocity profile for an inverse hybrid"
        " and returns the results of the velocity v as a function of the"
        " self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_mixed_hybrid_minus,
        v_plus,
        xi_wall,
        "Compute the velocity profile for a mixed hybrid with v- = c-.\n\n"
        "This function solves the fluid velocity profile for mixed hybrids"
        " with v- = c- and returns the results of the velocity v as a function"
        " of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    BIND_COMPUTE_VELOCITY_PROFILE(
        compute_velocity_profile_mixed_hybrid_plus,
        xi_wall,
        v_minus,
        "Compute the velocity profile for a mixed hybrid with v+ = c+.\n\n"
        "This function solves the fluid velocity profile for mixed hybrids"
        " with v+ = c+ and returns the results of the velocity v as a function"
        " of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "c_plus : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "c_minus : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    m.def(
        "get_shock_index",
        // throws exception if solution has no shock front -> no idx2obj needed
        get_shock_index, 
        py::arg("sol_type"), py::arg("size"),
        "Get the position of the shock front.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "size : int\n    Length of the profile vectors\n"
        "\nReturns\n-------\n"
        "idx : int\n    Index of the first element outside the shock front"
    );

    m.def(
        "get_shock_front",
        [](
            profile_calculator::HydroSolutionType sol_type, 
            py::array_t<double> xi, 
            py::array_t<double> v
        ) { 
            auto [xish, vsh] = get_shock_front(sol_type, np2vec(xi), np2vec(v));
            return py::make_tuple(xish, vsh);
        },
        py::arg("sol_type"), py::arg("xi"), py::arg("v"),
        "Get the coordinate and velocity of the shock front.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Self-similar coordinate xi=r/t of the profile\n"
        "v : numpy.ndarray, shape (N,)\n    Fluid velocity profile v(xi)\n"
        "\nReturns\n-------\n"
        "xi_sh : float\n    Self-similar coordinate of the shock front\n"
        "v_sh : float\n"
        "    Fluid velocity at the shock front (in the fluid rest frame)"
    );

    m.def(
        "compute_velocity_profile",
        [](
            double xi_wall, 
            double alpha_N, 
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
         ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile(
                xi, v, xi_wall, alpha_N, cp, cm, 
                recover_dynamic_default(step_size, profile_step_size), 
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("xi_wall"), 
        py::arg("alpha_N"), 
        py::arg("cp") = constants::cs_bag, 
        py::arg("cm") = constants::cs_bag,
        py::arg("step_size") = py::none(),
        py::arg("solve_mode") = py::none(),
        "Compute the velocity profile for a bubble in a cosmological phase"
        " transition.\n\n"
        "This function solves the fluid velocity profile for a cosmological"
        " phase transition and returns the results of the velocity v as a"
        " function of the self-similar coordinate xi.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "alpha_N : float\n    Pseudotrace phase transition strength alphaN"
        " (far in front of the wall)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile",
        [](
            HydroSolutionType sol_type, 
            double xi_wall, 
            double alpha_p, 
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile(
                sol_type, xi, v, xi_wall, alpha_p, cp, cm, 
                recover_dynamic_default(step_size, profile_step_size), 
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("sol_type"), 
        py::arg("xi_wall"),
        py::arg("alpha_p"), 
        py::arg("cp") = constants::cs_bag, 
        py::arg("cm") = constants::cs_bag,
        py::arg("step_size") = py::none(),
        py::arg("solve_mode") = py::none(),
        "Compute the velocity profile for a bubble in a cosmological phase"
        " transition.\n\n"
        "This function solves the fluid velocity profile for a cosmological"
        " phase transition and returns the results of the velocity v as a"
        " function of the self-similar coordinate xi.\n\n"
        "If a particular alpha+ root is already known, this overload"
        " constructs the corresponding stationary profile directly for the"
        " selected sol_type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "alpha_p : float\n    Pseudotrace strength alpha+ ahead of the wall\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : int\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    m.def(
        "integrate_enthalpy_segment",
        &integrate_enthalpy_segment,
        py::arg("xi1"),
        py::arg("xi2"),
        py::arg("v1"),
        py::arg("v2"),
        py::arg("cs") = constants::cs_bag,
        "Calculate a segment of the enthalpy density integral.\n\n"
        "This function calculates the segment of the integral in Eq. (2.12) of"
        " [Barni+ (2027)] assuming that xi is linear in v within the segment"
        " (xi1,xi2), i.e. xi = xi0 + v d(xi)/d(v), where"
        " xi0 = (xi1 v2 - xi2 v1)/(v2 - v1) and d(xi)/d(v) = (xi2-xi1)/(v2-v1)."
        " The integral can then be calculated using partial fraction"
        " decomposition.\n"
        "\nParameters\n----------\n"
        "xi1 : float\n"
        "    Self-similar coordinate xi1 at the beginning of the segment\n"
        "xi2 : float\n"
        "    Self-similar coordinate xi2 at the end of the segment\n"
        "v1 : float\n"
        "    Fluid velocity v1 at the beginning of the segment\n"
        "v2 : float\n"
        "    Fluid velocity v2 at the end of the segment\n"
        "cs : float, optional (default: cs_bag)\n    Speed of sound cs\n"
        "\nReturns\n-------\n"
        "I : float\n    Value of the integral I for the segment (xi1,xi2)"
    );

    m.def(
        "compute_enthalpy_density",
        [](
            py::array_t<double> xi, 
            py::array_t<double> v, 
            size_t idx_wall, 
            py::object idx_shock,
            double cp,
            double cm
        ) { 
            return vec2np(compute_enthalpy_density(
                np2vec(xi), np2vec(v), idx_wall, obj2idx(idx_shock), cp, cm
            ));
        },
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("idx_wall"), 
        py::arg("idx_shock")=py::none(),
        py::arg("cp") = constants::cs_bag, 
        py::arg("cm") = constants::cs_bag,
        "Compute the enthalpy density profile.\n\n"
        "Compute the enthalpy density profile from the fluid velocity profile,"
        " Eq. (2.12) of [Barni+ (2027)].\n"
        "\nParameters\n----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "idx_wall : int\n"
        "    Index of the first element outside the bubble\n"
        "idx_shock : int or None, optional (default: None)\n"
        "    Index of the first element outside the shock front (use None if"
        " there is no shock front)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Enthalpy density profile w(xi)"
    );

    m.def(
        "compute_energy_density_profile",
        [](
            py::array_t<double> w, 
            size_t idx_wall, 
            double epsilon,
            double cp,
            double cm
        ) { 
            return vec2np(compute_energy_density_profile(
                np2vec(w), idx_wall, epsilon, cp, cm
            )); 
        },
        py::arg("w"),
        py::arg("idx_wall"),
        py::arg("epsilon"),
        py::arg("cp") = constants::cs_bag, 
        py::arg("cm") = constants::cs_bag,
        "Compute the energy density profile.\n\n"
        "We assume the constant-sound-speed (or mu-nu) model,"
        " cf. [Leitao+ (2015)] and [Giese+ (2020)], i.e."
        " e(xi) = w(xi)/(1+cs^2(xi)) + epsilon(xi), where"
        " cs(xi) = c+ Theta(xi-xiw) + c- Theta(xiw-xi) and"
        " epsilon(xi) = epsilon+ Theta(xi-xiw) + epsilon- Theta(xiw-xi) is the"
        " vacuum energy contribution. Given the bag constant epsilon ="
        " epsilon+ - epsilon-, we set epsilon+ = epsilon, epsilon- = 0 if"
        " epsilon >= 0, and epsilon+ = 0, epsilon- = -epsilon otherwise.\n"
        "\nParameters\n----------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : int\n"
        "    Index of the first element outside the bubble\n"
        "epsilon : float\n    Bag constant epsilon = epsilon+ - epsilon-\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "e : numpy.ndarray, shape (N,)\n    Energy density profile e(xi)"
    );

    m.def(
        "compute_pressure_profile",
        [](
            py::array_t<double> w,
            size_t idx_wall,
            double epsilon,
            double cp,
            double cm
        ) { 
            return vec2np(compute_pressure_profile(
                np2vec(w), idx_wall, epsilon, cp, cm
            )); 
        },
        py::arg("w"),
        py::arg("idx_wall"),
        py::arg("epsilon"),
        py::arg("cp") = constants::cs_bag, 
        py::arg("cm") = constants::cs_bag,
        "Compute the pressure profile.\n\n"
        "We assume the constant-sound-speed (or mu-nu) model,"
        " cf. [Leitao+ (2015)] and [Giese+ (2020)], i.e."
        " p(xi) = w(xi)*cs^2(xi)/(1+cs^2(xi)) - epsilon(xi), where"
        " cs(xi) = c+ Theta(xi-xiw) + c- Theta(xiw-xi) and"
        " epsilon(xi) = epsilon+ Theta(xi-xiw) + epsilon- Theta(xiw-xi) is the"
        " vacuum energy contribution. Given the bag constant epsilon ="
        " epsilon+ - epsilon-, we set epsilon+ = epsilon, epsilon- = 0 if"
        " epsilon >= 0, and epsilon+ = 0, epsilon- = -epsilon otherwise.\n"
        "\nParameters\n----------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : int\n"
        "    Index of the first element outside the bubble\n"
        "epsilon : float\n    Bag constant epsilon = epsilon+ - epsilon-\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "p : numpy.ndarray, shape (N,)\n    Pressure profile p(xi)"
    );

    m.def(
        "compute_temperature_profile",
        [](
            py::array_t<double> w,
            size_t idx_wall,
            double T_N,
            double Psi_N,
            double cp,
            double cm
        ) { 
            return vec2np(compute_temperature_profile(
                np2vec(w), idx_wall, T_N, Psi_N, cp, cm
            )); 
        },
        py::arg("w"),
        py::arg("idx_wall"),
        py::arg("T_N"),
        py::arg("Psi_N"),
        py::arg("cp") = constants::cs_bag, 
        py::arg("cm") = constants::cs_bag,
        "Compute the temperature profile.\n\n"
        "We assume the constant-sound-speed (or mu-nu) model,"
        " cf. [Leitao+ (2015)] and [Giese+ (2020)], i.e."
        " w = (1+cs^2) a T^nu, where nu = 1 + 1/cs^2 and a is proportional"
        " to the effective number of relativistic degrees of freedom in the"
        " phase. We use the equal-temperature enthalpy ratio across the wall"
        " at the nucleation temperature, Psi_N = w-(T_N)/w+(T_N), to match"
        " temperatures at the bubble wall.\n"
        "\nParameters\n----------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : int\n"
        "    Index of the first element outside the bubble\n"
        "T_N : float\n    Nucleation temperature (far ahead of the wall)\n"
        "Psi_N : float\n    Enthalpy ratio across the wall at nucleation\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "T : numpy.ndarray, shape (N,)\n    Temperature profile T(xi)"
    );

    m.def(
        "compute_entropy_density_profile",
        [](
            py::array_t<double> w,
            py::array_t<double> T
        ) { 
            return vec2np(compute_entropy_density_profile(
                np2vec(w), np2vec(T)
            )); 
        },
        py::arg("w"),
        py::arg("T"),
        "Compute the entropy density profile.\n"
        "\nParameters\n----------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "T : numpy.ndarray, shape (N,)\n"
        "    Array containing the temperature T\n"
        "\nReturns\n-------\n"
        "s : numpy.ndarray, shape (N,)\n    Entropy density profile s(xi)"
    );

    m.def(
        "check_entropy_flux",
        [](
            py::array_t<double> xi,
            py::array_t<double> v,
            py::array_t<double> s
        ) -> py::object {
            return idx2obj(check_entropy_flux(
                np2vec(xi), np2vec(v), np2vec(s)
            ));
        },
        py::arg("xi"),
        py::arg("v"),
        py::arg("s"),
        "Check the entropy flux at discontinuities of a profile.\n\n"
        "Check that the entropy flux at a discontinuity does not violate the"
        " second law of thermodynamics, i.e. that it does not decrease in the"
        " rest frame of the discontinuity, s- gamma- v- >= s+ gamma+ v+.\n\n"
        "Note: The entropy density depends on the choice of Psi_N, cf."
        " compute_temperature_profile.\n"
        "\nParameters\n----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "s : numpy.ndarray, shape (N,)\n"
        "    Array containing the entropy density s\n"
        "\nReturns\n-------\n"
        "idx : int or None\n"
        "    Index of the first element outside the first discontinuity that"
        " violates the condition, or None if the profile is fine"
    );

    m.def(
        "compute_alpha_bag_from_pseudotrace",
        &compute_alpha_bag_from_pseudotrace,
        py::arg("alpha"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Compute the bag strength parameter from the pseudotrace strength"
        " parameter.\n"
        "\nParameters\n----------\n"
        "alpha : float\n    Pseudotrace transition strength\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "alpha_bag : float\n    Bag strength parameter"
    );

    m.def(
        "compute_alpha_pseudotrace_from_bag",
        &compute_alpha_pseudotrace_from_bag,
        py::arg("alpha"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Compute the pseudotrace strength parameter from the bag strength"
        " parameter.\n"
        "\nParameters\n----------\n"
        "alpha : float\n    Bag strength parameter\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "alpha : float\n    Pseudotrace transition strength"
    );

    m.def(
        "rescale_alpha_by_enthalpy",
        &rescale_alpha_by_enthalpy,
        py::arg("alpha"),
        py::arg("r_w"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Compute the pseudotrace strength parameter at some point (in front of"
        " the wall) from the value at some other point and the enthalpy ratio."
        "\n"
        "\nParameters\n----------\n"
        "alpha : float\n    Pseudotrace transition strength alpha_i\n"
        "r_w : float\n    Enthalpy ratio r_w = w_i/w_f\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "alpha : float\n    Pseudotrace transition strength alpha_f"
    );

    m.def(
        "compute_epsilon",
        &compute_epsilon,
        py::arg("w"), 
        py::arg("alpha"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        "Compute the bag constant (vacuum energy density).\n"
        "\nParameters\n----------\n"
        "w : float\n    Enthalpy density, either w+ (directly in front) or wN"
        " (far in front of the bubble wall)\n"
        "alpha : float\n    Pseudotrace transition strength, either alpha+"
        " (directly in front) or alphaN (far in front of the bubble wall)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "\nReturns\n-------\n"
        "epsilon : float\n    Bag constant epsilon = epsilon+ - epsilon-"
    );

    m.def(
        "compute_alpha_N",
        [](
            HydroSolutionType sol_type,
            double alpha_plus,
            double xi_wall,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ) {
            return compute_alpha_N(
                sol_type, alpha_plus, xi_wall, cp, cm,
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        },
        py::arg("sol_type"),
        py::arg("alpha_plus"),
        py::arg("xi_wall"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute the pseudotrace transition strength alphaN from the value in"
        " front of the bubble wall alpha+.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "alpha_plus : float\n    Pseudotrace transition strength alpha+"
        " immediately in front of the bubble wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "alpha_N : float\n"
        "    Pseudotrace transition strength alphaN far in front of the"
        " bubble wall"
    );

    m.def(
        "find_alpha_plus",
        [](
            double alpha_N,
            double xi_wall,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ) {
            return find_alpha_plus(
                alpha_N, xi_wall, cp, cm, 
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        },
        py::arg("alpha_N"),
        py::arg("xi_wall"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Find the pseudotrace transition strength alpha+ from the value alphaN"
        " in the false vacuum phase.\n\n"
        "The solution type is determined by classify_wave_fluid_frame. The"
        " automatic classification may report a non-injective sound-speed"
        " region depending on the noninjectivity check mode (see"
        " set_noninjectivity_check).\n\n"
        "This is a single-root inversion: it returns one alpha+. On a folded"
        " branch, more than one alpha+ can correspond to the same alphaN even"
        " within the same hydrodynamic branch.\n"
        "\nParameters\n----------\n"
        "alpha_N : float\n"
        "    Pseudotrace transition strength alphaN far in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "alpha_plus : float\n    Value of the pseudotrace strength alpha+"
    );

    m.def(
        "find_alpha_plus",
        [](
            HydroSolutionType sol_type,
            double alpha_N,
            double xi_wall,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ) {
            return find_alpha_plus(
                sol_type, alpha_N, xi_wall, cp, cm, 
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        },
        py::arg("sol_type"),
        py::arg("alpha_N"),
        py::arg("xi_wall"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Find the pseudotrace transition strength alpha+ from the value alphaN"
        " in the false vacuum phase.\n\n"
        "The search is restricted to the selected sol_type. This is a"
        " single-root inversion: it returns one alpha+. On a folded branch,"
        " more than one alpha+ can correspond to the same alphaN even within"
        " the same sol_type; selecting a branch explicitly does not make this"
        " function an all-roots solver.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "alpha_N : float\n"
        "    Pseudotrace transition strength alphaN far in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "alpha_plus : float\n    Value of the pseudotrace strength alpha+"
    );

    m.def(
        "compute_efficiency",
        [](
            py::array_t<double> xi, 
            py::array_t<double> v, 
            py::array_t<double> w, 
            size_t idx_wall, 
            double alpha_N,
            double cp,
            double cm,
            EfficiencyNormalization normalization
        ) {
            return compute_efficiency(
                np2vec(xi),
                np2vec(v),
                np2vec(w),
                idx_wall,
                alpha_N,
                cp,
                cm,
                normalization
            );
        },
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("w"), 
        py::arg("idx_wall"), 
        py::arg("alpha_N"),
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("normalization")=EfficiencyNormalization::TotalEnergyDensity,
        "Compute the efficiency.\n\n"
        "Compute the kinetic energy fraction K = rho_kin/rho_tot or the"
        " efficiency factor kappa_sw for converting vacuum energy to bulk"
        " kinetic energy, Eq. (6.3) and (6.4) of [Barni+ (2027)]. See also"
        " Eq. (2.30) of [Espinosa+ (2010)] or Eq. (4.10) of [Barni+ (2024)].\n"
        "\nParameters\n----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : int\n"
        "    Index of the first element outside the bubble\n"
        "alpha_N : float\n"
        "    Pseudotrace transition strength alphaN (far in front of the"
        " wall)\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "normalization : EfficiencyNormalization, optional (default:"
        " EfficiencyNormalization.TotalEnergyDensity)\n"
        "    EfficiencyNormalization switch to determine how to normalize the"
        " efficiency\n"
        "\nReturns\n-------\n"
        "kappa_or_K : float\n"
        "    Efficiency kappa_sw or K (depending on normalization)"
    );

    m.def(
        "compute_enthalpy_ratio_LTE",
        [](
            double xi_wall,
            double alpha_N,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ) {
            return compute_enthalpy_ratio_LTE(
                xi_wall, alpha_N, cp, cm, 
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        },
        py::arg("xi_wall"), 
        py::arg("alpha_N"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Compute the enthalpy ratio Psi+ in local thermal equilibrium.\n\n"
        "Compute the enthalpy ratio Psi+ = w-(T+)/w+(T+) ahead of and behind"
        " the wall, evaluated at the temperature ahead of the wall. This"
        " assumes local thermal equilibrium and uses entropy"
        " conservation to relate the temperatures on both sides of the wall,"
        " i.e. s+ v+ gamma+ = s- v- gamma-.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Pseudotrace transition strength alphaN far in front of the wall\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "Psi : float\n    Enthalpy ratio Psi+ across the wall"
    );

    m.def(
        "find_LTE_wall_velocity",
        [](
            double alpha_N,
            double Psi,
            double xi_min,
            double xi_max,
            double cp,
            double cm,
            py::object step_size,
            py::object solve_mode
        ){
            return find_LTE_wall_velocity(
                alpha_N, Psi, xi_min, xi_max, cp, cm,
                recover_dynamic_default(step_size, profile_step_size),
                recover_dynamic_default(solve_mode, profile_solve_mode)
            );
        },
        py::arg("alpha_N"), 
        py::arg("Psi"), 
        py::arg("xi_min"), 
        py::arg("xi_max"), 
        py::arg("cp")=constants::cs_bag,
        py::arg("cm")=constants::cs_bag,
        py::arg("step_size")=py::none(),
        py::arg("solve_mode")=py::none(),
        "Find the LTE wall velocity xi_w for a given enthalpy ratio Psi+.\n\n"
        "Find the wall velocity xi_w for a given enthalpy ratio Psi+ = w-(T+)/"
        "w+(T+) assuming local thermal equilibrium.\n"
        "\nParameters\n----------\n"
        "alpha_N : float\n"
        "    Pseudotrace transition strength alphaN far in front of the wall\n"
        "Psi : float\n    Enthalpy ratio Psi+\n"
        "xi_min : float\n    Minimal wall velocity xi_w_min to consider\n"
        "xi_max : float\n    Maximal wall velocity xi_w_max to consider\n"
        "cp : float, optional (default: cs_bag)\n"
        "    Speed of sound c+ ahead of the wall\n"
        "cm : float, optional (default: cs_bag)\n"
        "    Speed of sound c- behind the wall\n"
        "step_size : float, optional (default: profile_calculator.profile_step_"
        "size)\n    Step size for the fluid velocity (or coordinate) when"
        " solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: profile_calculator."
        "profile_solve_mode)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "xi_wall : float\n    Wall velocity xi_w"
    );

    m.def(
        "check_profiles",
        [](
            py::array_t<double> xi, 
            py::array_t<double> profile, 
            bool throw_error
        ) {
            return check_profiles(np2vec(xi), np2vec(profile), throw_error);
        },
        py::arg("xi"), py::arg("profile"), py::arg("throw_error")=false,
        "Check profile for consistency.\n\n"
        "This function returns True if xi and profile have the same size and at"
        " least size 6. Otherwise it throws an error (if throw_error is True)"
        " or returns False.\n"
        "\nParameters\n----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Self-similar coordinate xi=r/t of the profile\n"
        "profile : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid profile\n"
        "throw_error : bool, optional (default: False)\n"
        "    Whether to throw an error if the test fails\n"
        "\nReturns\n-------\n"
        "valid : bool\n    Whether the profile passes the test"
    );

    // getter and setter for settings:

    m.def(
        "get_equal_rel_tol", 
        []() { return equal_rel_tol; }, 
        "Get the relative tolerance for when to consider velocities/coordinates"
        " as equal."
    );
    m.def(
        "set_equal_rel_tol", 
        [](double rel_tol) { equal_rel_tol = rel_tol; }, 
        py::arg("rel_tol"), 
        "Set the relative tolerance for when to consider velocities/coordinates"
        " as equal."
    );

    m.def(
        "get_equal_abs_tol", 
        []() { return equal_abs_tol; }, 
        "Get the absolute tolerance for when to consider velocities/coordinates"
        " as equal."
    );
    m.def(
        "set_equal_abs_tol", 
        [](double abs_tol) { equal_abs_tol = abs_tol; }, 
        py::arg("abs_tol"), 
        "Set the absolute tolerance for when to consider velocities/coordinates"
        " as equal."
    );

    m.def(
        "get_profile_step_size",
        []() { return profile_step_size; },
        "Get the default step size when integrating the fluid profile."
    );
    m.def(
        "set_profile_step_size",
        [](double step_size) { profile_step_size = step_size; },
        py::arg("step_size"),
        "Set the default step size when integrating the fluid profile."
    );

    m.def(
        "get_profile_solve_mode",
        []() { return profile_solve_mode; },
        "Get the default profile integration method."
    );
    m.def(
        "set_profile_solve_mode",
        [](ProfileSolveMode solve_mode) { profile_solve_mode = solve_mode; },
        py::arg("solve_mode"),
        "Set the default profile integration method."
    );

    m.def(
        "get_alpha_rel_tol", 
        []() { return alpha_rel_tol; }, 
        "Get the relative tolerance when finding alpha_+ from alpha_N."
    );
    m.def(
        "set_alpha_rel_tol", 
        [](double rel_tol) { alpha_rel_tol = rel_tol; }, 
        py::arg("rel_tol"), 
        "Set the relative tolerance when finding alpha_+ from alpha_N."
    );

    m.def(
        "get_alpha_abs_tol", 
        []() { return alpha_abs_tol; }, 
        "Get the absolute tolerance when finding alpha_+ from alpha_N."
    );
    m.def(
        "set_alpha_abs_tol", 
        [](double abs_tol) { alpha_abs_tol = abs_tol; }, 
        py::arg("abs_tol"), 
        "Set the absolute tolerance when finding alpha_+ from alpha_N."
    );

    m.def(
        "get_xi_wall_rel_tol", 
        []() { return xi_wall_rel_tol; }, 
        "Get the relative tolerance when finding the wall velocity from the"
        " LTE condition."
    );
    m.def(
        "set_xi_wall_rel_tol", 
        [](double rel_tol) { xi_wall_rel_tol = rel_tol; }, 
        py::arg("rel_tol"), 
        "Set the relative tolerance when finding the wall velocity from the"
        " LTE condition."
    );

    m.def(
        "get_xi_wall_abs_tol", 
        []() { return xi_wall_abs_tol; }, 
        "Get the absolute tolerance when finding the wall velocity from the"
        " LTE condition."
    );
    m.def(
        "set_xi_wall_abs_tol", 
        [](double abs_tol) { xi_wall_abs_tol = abs_tol; }, 
        py::arg("abs_tol"), 
        "Set the absolute tolerance when finding the wall velocity from the"
        " LTE condition."
    );
}