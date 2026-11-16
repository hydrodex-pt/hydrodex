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

#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "hydrodex/constants.hpp"

/**
 * @namespace hydrodex::profile_calculator
 * @brief Functions to calculate the fluid profile of a bubble in a cosmological
 *     phase transition.
 * 
 * This is based on 
 * 1. @anchor ref_Barni_2027
 *    G. Barni and E. Madge,
 *    "HydróDex: gotta match 'em all!"
 *    [arXiv:2611.17411 [physics.flu-dyn]](https://arxiv.org/abs/2611.17411).
 * 2. @anchor ref_Barni_2026
 *    G. Barni, S. Blasi, E. Madge and M. Vanvlasselaer,
 *    "Gravitational waves from the sound shell model: direct and inverse phase
 *    transitions in the early Universe,"
 *    [JCAP **3** (2026), 066](https://doi.org/10.1088/1475-7516/2026/03/066),
 *    [arXiv:2510.21439 [hep-ph]](https://arxiv.org/abs/2510.21439).
 * 3. @anchor ref_Barni_2025
 *    G. Barni, S. Blasi and M. Vanvlasselaer, 
 *    "Inverse bubbles from broken symmetries,"
 *    [PRD **112** (2025), 9](https://doi.org/10.1103/pmf7-vlvc),
 *    [arXiv:2503.01951 [hep-ph]](https://arxiv.org/abs/2503.01951).
 * 4. @anchor ref_Barni_2024
 *    G. Barni, S. Blasi and M. Vanvlasselaer, 
 *    "The hydrodynamics of inverse phase transitions,"
 *    [JCAP **10** (2024), 042](https://doi.org/10.1088/1475-7516/2024/10/042),
 *    [arXiv:2406.01596 [hep-ph]](https://arxiv.org/abs/2406.01596).
 * 5. @anchor ref_Giese_2020
 *    F. Giese, T. Konstandin and J. van de Vis,
 *    "Model-independent energy budget of cosmological first-order phase
 *    transitions - A sound argument to go beyond the bag model,"
 *    [JCAP **7** (2020), 057](https://doi.org/10.1088/1475-7516/2020/07/057),
 *    [arXiv:2004.06995 [astro-ph.CO]](https://arxiv.org/abs/2004.06995).
 * 6. @anchor ref_Leitao_2015
 *    L. Leitao and A. Megevand,
 *    "Hydrodynamics of phase transition fronts and the speed of sound in the
 *    plasma,"
 *    [NPB **891** (2015), 159-199]
 *    (https://doi.org/10.1016/j.nuclphysb.2014.12.008),
 *    [arXiv:1410.3875 [hep-ph]](https://arxiv.org/abs/1410.3875).
 * 7. @anchor ref_Espinosa_2010
 *    J. R. Espinosa and J. M. No,
 *    "Energy Budget of Cosmological First-order Phase Transitions,"
 *    [JCAP **06** (2010), 028](https://doi.org/10.1088/1475-7516/2010/06/028),
 *    [arXiv:1004.4187 [hep-ph]](https://arxiv.org/abs/1004.4187).
 */
namespace hydrodex::profile_calculator {

/// @brief Types of the hydrodynamic solutions.
enum class HydroSolutionType { 
    /// weak (subsonic) direct deflagration
    Deflagration,
    /// weak (supersonic) direct detonation
    Detonation,
    /// direct hybrid solution (supersonic deflagration) with \f$v_+\f$ on '(-)'
    /// branch
    Hybrid,
    /// direct hybrid solution (supersonic deflagration) with \f$v_+\f$ on '(+)'
    /// branch
    HybridPlus,
    /// weak (supersonic) inverse deflagration
    InverseDeflagration, 
    /// weak (subsonic) inverse detonation with \f$v_+\f$ on '(-)' branch
    InverseDetonation,   
    /// weak (subsonic) inverse detonation with \f$v_+\f$ on '(+)' branch
    InverseDetonationPlus,   
    /// inverse hybrid solution (subsonic inverse deflagration)
    InverseHybrid,
    /// mixed direct-inverse hybrid solution with \f$v_-=c_-\f$
    MixedHybridMinus,
    /// mixed direct-inverse hybrid solution with \f$v_+=c_+\f$
    MixedHybridPlus
};

/// @brief Normalization convention for the hydrodynamic efficiency.
enum class EfficiencyNormalization {
    /**
     * @brief Normalize to the bubble volume.
     * 
     * This corresponds to the kinetic energy density \f$\rho_\mathrm{kin}\f$.
     */
    BubbleVolume,
    /**
     * @brief Normalize to the total energy density.
     * 
     * This corresponds to the kinetic energy fraction
     * \f$K = \frac{\rho_\mathrm{kin}}{\rho_\mathrm{tot}}\f$.
     */
    TotalEnergyDensity,
    /**
     * @brief Normalize to the vacuum energy density.
     * 
     * This corresponds to
     * \f$\kappa = \frac{\rho_\mathrm{kin}}{\rho_\mathrm{vac}}\f$
     * with \f$\rho_\mathrm{vac} = \frac{4}{3}\, \pi\,\xi_w^3\,\epsilon\f$.
     * 
     * For negative \f$\epsilon < 0\f$, this option normalizes to the initial
     * thermal energy \f$\rho_\mathrm{ini}\f$ instead, where
     * \f$\rho_\mathrm{ini} = \frac{4}{3}\,\pi\,\xi_w^3\,\frac{w_N}{1+c_+^2}\f$.
     */
    VacuumEnergy,
    /**
     * @brief Normalize to the pseudotrace strength parameter
     *     \f$\alpha_N\f$ at nucleation.
     * 
     * This corresponds to \f$\bar\kappa =
     * \frac{4\,\rho_\mathrm{kin}}{3\,\alpha_N w_N}\f$.
     * 
     * For negative \f$\alpha_N < 0\f$, this option normalizes to
     * the initial thermal energy \f$\rho_\mathrm{ini}\f$ instead, where
     * \f$\rho_\mathrm{ini} = \frac{4}{3}\,\pi\,\xi_w^3\,\frac{w_N}{1+c_+^2}\f$.
     */
    AlphaPseudoTrace
};

/** 
 * @brief Selection which differential equation to solve when calculating the
 *     velocity profile.
 */
enum class ProfileSolveMode {
    D_XI_D_V,           ///< always solve \f$\frac{d\xi}{dv}\f$
    D_V_D_XI,           ///< always solve \f$\frac{dv}{d\xi}\f$
    /** 
     * solve \f$\frac{dv}{d\xi}\f$ when the solution is expected to end in a
     * shock and \f$\frac{d\xi}{dv}\f$ otherwise 
     */
    D_V_D_XI_FOR_SHOCKS 
};

/**
 * @brief Branch selector for the sign in front of the square root in the
 *     matching for \f$v_\pm\f$.
 */
enum class RootBranch {
    Plus, ///< branch with '+' sign in front of the square root
    Minus, ///< branch with '-' sign in front of the square root
};

/// relative tolerance when considering velocities/coordinates as equal
extern double equal_rel_tol;
/// absolute tolerance when considering velocities/coordinates as equal
extern double equal_abs_tol;
/// default step size when integrating the fluid profile
extern double profile_step_size;
/// default profile integration method
extern ProfileSolveMode profile_solve_mode;
/// relative tolerance when finding alpha_+ from alpha_N
extern double alpha_rel_tol; 
/// absolute tolerance when finding alpha_+ from alpha_N
extern double alpha_abs_tol; 
/// relative tolerance when finding wall velocity from LTE condition
extern double xi_wall_rel_tol;
/// absolute tolerance when finding wall velocity from LTE condition 
extern double xi_wall_abs_tol;

/**
 * @brief Convert \c HydroSolutionType \c enum value to a human-readable string.
 * @param type The \c HydroSolutionType value
 * @return String representation of the enum value
 */
const std::string to_string(HydroSolutionType type);

/**
 * @brief Stream output operator for \c HydroSolutionType.
 * @param os The output stream
 * @param type The \c HydroSolutionType value to print
 * @return The output stream with the string representation of the
 *     \c HydroSolutionType appended
 */
std::ostream& operator<<(std::ostream& os, HydroSolutionType type);

/**
 * @brief Left multiplication operator for \c RootBranch.
 * 
 * Multiplies a number by the sign corresponding to the \c RootBranch value.
 * 
 * @param branch \c RootBranch value determining the sign
 * @param x Number to be multiplied
 * @return \p x for \c RootBranch::Plus and \f$-x\f$ for \c RootBranch::Minus
 */
template <typename T>
constexpr T operator*(RootBranch branch, T x) {
    return (branch == RootBranch::Plus) ? x : -x;
}

/**
 * @brief Right multiplication operator for \c RootBranch.
 * 
 * Multiplies a number by the sign corresponding to the \c RootBranch value.
 * 
 * @param x Number to be multiplied
 * @param branch \c RootBranch value determining the sign
 * @return \p x for \c RootBranch::Plus and \f$-x\f$ for \c RootBranch::Minus
 */
template <typename T>
constexpr T operator*(T x, RootBranch branch) {
    return branch * x;
}

/**
 * @brief Check if two numbers are equal.
 * 
 * This function checks whether the two numbers \p x and \p y are equal within
 * the relative precision \p eps_rel and the absolute precision \p eps_abs.
 * 
 * @param x First number
 * @param y Second number
 * @param eps_rel Relative precision in the comparison
 *     (default: \c profile_calculator::equal_rel_tol)
 * @param eps_abs Absolute precision in the comparison 
 *     (default: \c profile_calculator::equal_abs_tol)
 * @return Whether the two numbers are equal within the tolerance
 */
bool are_equal(
    double x, 
    double y, 
    double eps_rel=equal_rel_tol, 
    double eps_abs=equal_abs_tol
);

/**
 * @brief Helper function to shift a value within a range.
 * 
 * This is a helper function to avoid calculating a function directly at a
 * minimal or maximal allowed value \p x, the value of \p x is shifted by a 
 * small amount \p eps_abs towards \p y, making sure that the result is not
 * shifted by more than a fraction \p eps_rel of the distance between \p x 
 * and \p y.
 * 
 * This function is used internally when calculating the minimal and maximal
 * values of \f$\alpha_N\f$ from \f$\alpha_+\f$.
 * 
 * @param x First number (to be shifted)
 * @param y Second number
 * @param eps_rel Maximal allowed relative shift (default: 
 *     \c profile_calculator::equal_rel_tol)
 * @param eps_abs Absolute shift (default: 
 *     \c profile_calculator::equal_abs_tol)
 * @return Shifted value of \p x
 */
double shift_safe(
    double x, 
    double y, 
    double eps_rel=equal_rel_tol,
    double eps_abs=equal_abs_tol
);

/**
 * @brief Compute the Lorentz boost factor \f$\gamma=\frac{1}{\sqrt{1-v^2}}\f$.
 * @param v Velocity \f$v\f$
 * @return Boost factor \f$\gamma\f$
 */
double gamma(double v);

/**
 * @brief Lorentz boosted velocity.
 * 
 * This boosts the velocity \f$v\f$ to the frame moving with velocity \f$\xi\f$,
 * Eq. (2.7) of @ref ref_Barni_2027 "Barni and Madge (2027)".
 * \f[ \mu(\xi,v) = \frac{\xi-v}{1-\xi v} \f]
 * 
 * @param xi Relative velocity \f$\xi\f$ of the new frame to which the velocity
 *     is boosted
 * @param v Velocity \f$v\f$ in the original frame
 * @return Boosted velocity \f$\mu(\xi, v)\f$
 */
double mu(double xi, double v);

/**
 * @brief Check if the velocity is between 0 and 1.
 * @param v Velocity to check
 * @return @c true if the velocity lies between 0 and 1, @c false otherwise
 */
bool check_velocity(double v);

/**
 * @brief Get the temperature exponent \f$\nu\f$ from the speed of sound
 *     \f$c_s\f$ in the constant-sound-speed model.
 * @param cs Speed of sound \f$c_s\f$
 * @return Temperature exponent \f$\nu\f$
 */
double get_nu_from_cs(double cs);

/**
 * @brief Get the speed of sound \f$c_s\f$ from the temperature exponent
 *     \f$\nu\f$ in the constant-sound-speed model.
 * @param nu Temperature exponent \f$\nu\f$
 * @return Speed of sound \f$c_s\f$
 */
double get_cs_from_nu(double nu);

/**
 * @brief Calculate the distance to the shock front.
 * 
 * This function returns the shock condition
 * \f$\mu(\xi_\text{sh}, v_\text{sh}) \xi_\text{sh} - c_s^2\f$.
 * It returns zero if \f$(\xi,v)\f$ corresponds to a shock front, a negative
 * value if the point is ahead of a shock front and a positive value if it is
 * behind the shock front.
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$
 * @param v Fluid velocity \f$v\f$
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return Distance to the shock front
 */
double shock_condition(double xi, double v, double cs=constants::cs_bag);

/**
 * @brief Calculate the distance to the Jouguet front.
 * 
 * This function returns the Jouguet condition
 * \f$\mu(\xi_\text{J}, v_\text{J})^2 - c_s^2\f$.
 * It returns zero if \f$(\xi,v)\f$ corresponds to a Jouguet front, a negative
 * value if the point is ahead of a Jouguet front and a positive value if it is
 * behind the Jouguet front.
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$
 * @param v Fluid velocity \f$v\f$
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return Distance to the Jouguet front
 */
double jouguet_condition(double xi, double v, double cs=constants::cs_bag);

/** 
 * @brief Calculate \f$v_+\f$ from \f$v_-\f$ and \f$\alpha_+\f$.
 * 
 * This function calculates the fluid velocity \f$v_+\f$ ahead of the bubble (or
 * shock) front from the velocity \f$v_-\f$ behind the front and the
 * pseudotrace transition strength parameter \f$\alpha_+\f$ ahead of the front. 
 * 
 * This function returns both branches of the solution for \f$v_+\f$.
 * 
 * @param vm Fluid velocity \f$v_-\f$ behind the front
 * @param alpha Pseudotrace transition strength parameter \f$\alpha_+\f$ ahead
 *     of the front
 * @param cm Speed of sound \f$c_-\f$ behind the front (default:
 *     \c constants::cs_bag)
 * @return The solutions on the two branches for the fluid velocity \f$v_+\f$
 *     ahead of the front
 */
std::pair<double, double> get_v_plus_from_v_minus(
    double vm, 
    double alpha,
    double cm = constants::cs_bag
);

/** 
 * @brief Calculate \f$v_+\f$ from \f$v_-\f$ and \f$\alpha_+\f$.
 * 
 * This function calculates the fluid velocity \f$v_+\f$ ahead of the bubble (or
 * shock) front from the velocity \f$v_-\f$ behind the front and the
 * pseudotrace transition strength parameter \f$\alpha_+\f$ ahead of the front. 
 * 
 * This function returns the '(-)' or '(+)' branch of the solution for \f$v_+\f$
 * depending on the value of \p branch.
 * 
 * @param vm Fluid velocity \f$v_-\f$ behind the front
 * @param alpha Pseudotrace transition strength parameter \f$\alpha_+\f$ ahead
 *     of the front
 * @param branch \c RootBranch selector to determine which branch to return
 * @param cm Speed of sound \f$c_-\f$ behind the front (default:
 *     \c constants::cs_bag)
 * @return Solution for the fluid velocity \f$v_+\f$ ahead of the front on the
 *     selected branch
 */
double get_v_plus_from_v_minus(
    double vm, 
    double alpha,
    RootBranch branch,
    double cm = constants::cs_bag
);

/** 
 * @brief Calculate \f$v_-\f$ from \f$v_+\f$ and \f$\alpha_+\f$.
 * 
 * This function calculates the fluid velocity \f$v_-\f$ behind the bubble (or
 * shock) front from the velocity \f$v_+\f$ and the pseudotrace transition
 * strength \f$\alpha_+\f$ ahead of the front. 
 * 
 * This function returns both branches of the solution for \f$v_-\f$.
 * 
 * @param vp Fluid velocity \f$v_+\f$ ahead of the front
 * @param alpha Pseudotrace transition strength parameter
 *     \f$\alpha_+\f$ ahead of the front
 * @param cm Speed of sound \f$c_-\f$ behind the front (default:
 *     \c constants::cs_bag)
 * @return The solutions on the two branches for the fluid velocity \f$v_-\f$
 *     behind the front
 */
std::pair<double, double> get_v_minus_from_v_plus(
    double vp,
    double alpha,
    double cm = constants::cs_bag
);

/** 
 * @brief Calculate \f$v_-\f$ from \f$v_+\f$ and \f$\alpha_+\f$.
 * 
 * This function calculates the fluid velocity \f$v_-\f$ behind the bubble (or
 * shock) front from the velocity \f$v_+\f$ and the pseudotrace transition
 * strength \f$\alpha_+\f$ ahead of the front. 
 * 
 * This function returns the '(-)' or '(+)' branch of the solution for \f$v_-\f$
 * depending on the value of \p branch.
 * 
 * @param vp Fluid velocity \f$v_+\f$ ahead of the front
 * @param alpha Pseudotrace transition strength parameter \f$\alpha_+\f$ ahead
 *     of the front
 * @param branch \c RootBranch selector to determine which branch to return
 * @param cm Speed of sound \f$c_-\f$ behind the front (default:
 *     \c constants::cs_bag)
 * @return Solution for the fluid velocity \f$v_-\f$ behind the front on the
 *     selected branch
 */
double get_v_minus_from_v_plus(
    double vp,
    double alpha,
    RootBranch branch,
    double cm = constants::cs_bag
);

/**
 * @brief Compute \f$\alpha_+\f$ from \f$v_+\f$ and \f$v_-\f$.
 * 
 * This function calculates the pseudotrace transition strength ahead of the
 * wall, \f$\alpha_+\f$, from the fluid velocities \f$v_+\f$ and \f$v_-\f$ in
 * the wall rest frame.
 * 
 * @note For \f$c_- = c_{s,\text{bag}} = \frac{1}{\sqrt{3}}\f$, this coincides
 *     with the bag-model transition strength \f$\alpha_{\text{bag}}^+\f$.
 * 
 * @param vp Fluid velocity \f$v_+\f$ ahead of the wall
 * @param vm Fluid velocity \f$v_-\f$ behind the wall
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Pseudotrace phase transition strength \f$\alpha_+\f$ ahead of the
 *     bubble wall
 */
double compute_alpha_plus(double vp, double vm, double cm=constants::cs_bag);

/**
 * @brief Get the extremum position \f$v_{+,\text{extr}}\f$.
 * 
 * The '(-)' and '(+)' branch of \f$v_-^{(\pm)}(v_+,\alpha_+)\f$ exhibit a
 * maximum or minimum, respectively, for strength parameters
 * \f$-\frac{1-c_-^2}{6\,c_-^2} < \alpha_+ < 0\f$.
 * 
 * This function returns the position \f$v_{+,\text{extr}}\f$ of the extremum,
 * cf. Eq. (3.9a) of @ref ref_Barni_2027 "Barni and Madge (2027)", with the
 * value of \f$\alpha_+\f$ such that the maximal or minimal value is given by
 * the input parameter \p vm, i.e.
 * \f$v_-^{(\pm)}(v_{+,\text{extr}},\alpha_+) = v_-\f$.
 * 
 * @param vm Value \f$v_{-,\text{extr}}\f$ of \f$v_-\f$ at the extremum
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Position \f$v_{+,\text{extr}}\f$ of the extremum
 */
double get_v_plus_extremum_from_v_minus(double vm, double cm=constants::cs_bag);

/**
 * @brief Get the extremum position \f$v_{+,\text{extr}}\f$.
 * 
 * The '(-)' and '(+)' branch of \f$v_-^{(\pm)}(v_+,\alpha_+)\f$ exhibit a
 * maximum or minimum, respectively, for strength parameters
 * \f$-\frac{1-c_-^2}{6\,c_-^2} < \alpha_+ < 0\f$.
 * 
 * This function returns the position \f$v_{+,\text{extr}}\f$ of the extremum,
 * cf. Eq. (3.9a) of @ref ref_Barni_2027 "Barni and Madge (2027)", for a given
 * value of \f$\alpha_+\f$.
 * 
 * @param alpha Pseudotrace transition strength parameter
 *     \f$\alpha_+\f$ ahead of the front
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Position \f$v_{+,\text{extr}}\f$ of the extremum
 */
double get_v_plus_extremum_from_alpha(
    double alpha, 
    double cm=constants::cs_bag
);

/**
 * @brief Get the extremum value \f$v_{-,\text{extr}}\f$.
 * 
 * The '(-)' and '(+)' branch of \f$v_-^{(\pm)}(v_+,\alpha_+)\f$ exhibit a
 * maximum or minimum, respectively, for strength parameters
 * \f$-\frac{1-c_-^2}{6\,c_-^2} < \alpha_+ < 0\f$.
 * 
 * This function returns the value at the extremum
 * \f$v_{-,\text{extr}}^{(\pm)} = v_-^{(\pm)}(v_{+,\text{extr}}, \alpha_+)\f$
 * on the branch determined by \p branch, cf. Eq. (3.9b) of
 * @ref ref_Barni_2027 "Barni and Madge (2027)", with \f$\alpha_+\f$ such that
 * the extremum occurs at the position \f$v_{+,\text{extr}}\f$ given by the
 * input parameter \p vp.
 * 
 * @param vp Position \f$v_{+,\text{extr}}\f$ of the extremum
 * @param branch \c RootBranch selector to determine which branch to return
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Value \f$v^{(\pm)}_{-,\text{extr}}\f$ at the extremum
 */
double get_v_minus_extremum_from_v_plus(
    double vp,
    RootBranch branch,
    double cm=constants::cs_bag
);  

/**
 * @brief Get the extremum value \f$v_{-,\text{extr}}\f$.
 * 
 * The '(-)' and '(+)' branch of \f$v_-^{(\pm)}(v_+,\alpha_+)\f$ exhibit a
 * maximum or minimum, respectively, for strength parameters
 * \f$-\frac{1-c_-^2}{6\,c_-^2} < \alpha_+ < 0\f$.
 * 
 * This function returns the value at the extremum
 * \f$v_{-,\text{extr}}^{(\pm)} = v_-^{(\pm)}(v_{+,\text{extr}}, \alpha_+)\f$
 * on the branch determined by \p branch, cf. Eq. (3.9b) of
 * @ref ref_Barni_2027 "Barni and Madge (2027)", for a given value of
 * \f$\alpha_+\f$.
 * 
 * @param alpha Pseudotrace transition strength parameter \f$\alpha_+\f$ ahead
 *     of the front
 * @param branch \c RootBranch selector to determine which branch to return
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Value \f$v^{(\pm)}_{-,\text{extr}}\f$ at the extremum
 */
double get_v_minus_extremum_from_alpha(
    double alpha, 
    RootBranch branch,
    double cm=constants::cs_bag
);  

/**
 * @brief Compute the Jouguet velocity for a direct phase transition.
 * @param alpha Pseudotrace strength parameter \f$\alpha_+\f$
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Jouguet velocity \f$v_J^\text{direct}\f$
 */
double direct_jouguet_velocity(double alpha, double cm=constants::cs_bag);

/**
 * @brief Compute the Jouguet velocity for an inverse phase transition.
 * @param alpha Pseudotrace strength parameter \f$\alpha_+\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Inverse Jouguet velocity \f$v_J^\text{inverse}\f$
 */
double inverse_jouguet_velocity(
    double alpha, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Minimal wall velocity \f$\xi_w\f$ for a given solution type.
 * @param sol_type Type of the solution
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return The minimal wall velocity for the solution type @p sol_type
 */
double xi_min(
    HydroSolutionType sol_type,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag
);

/**
 * @brief Maximal wall velocity \f$\xi_w\f$ for a given solution type.
 * @param sol_type Type of the solution
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return The maximal wall velocity for the solution type @p sol_type
 */
double xi_max(
    HydroSolutionType sol_type,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag
);

/**
 * @brief Check if @p xi_w is a valid wall velocity for solutions of type
 *     @p sol_type.
 * @param sol_type Type of the solution
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @throws std::invalid_argument if the wall velocity is not valid
 */
void check_wall_velocity(
    HydroSolutionType sol_type,
    double xi_w,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     detonation solutions.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct detonations.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for direct detonations
 */
double alpha_max_detonation(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     detonation solutions.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct detonations.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for direct detonations
 */
double alpha_min_detonation(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     hybrid solutions on the '(+)' branch.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct hybrids on the '(+)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for direct hybrids (on
 *     the '(+)' branch)
 */
double alpha_max_hybrid_plus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     hybrid solutions on the '(+)' branch.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct hybrids on the '(+)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for direct hybrids (on
 *     the '(+)' branch)
 */
double alpha_min_hybrid_plus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     hybrid solutions on the '(-)' branch.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct hybrids on the '(-)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for direct hybrids (on
 *     the '(-)' branch)
 */
double alpha_max_hybrid(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     hybrid solutions on the '(-)' branch.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct hybrids on the '(-)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for direct hybrids (on
 *     the '(-)' branch)
 */
double alpha_min_hybrid(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     deflagration solutions.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct deflagrations.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for direct deflagrations
 */
double alpha_max_deflagration(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     deflagration solutions.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct deflagrations.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for direct deflagrations
 */
double alpha_min_deflagration(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse detonation solutions on the '(-)' branch.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse detonations on the '(-)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for inverse detonations
 *     (on the '(-)' branch)
 */
double alpha_max_inverse_detonation(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse detonation solutions on the '(-)' branch.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse detonations on the '(-)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for inverse detonations
 *     (on the '(-)' branch)
 */
double alpha_min_inverse_detonation(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse detonation solutions on the '(+)' branch.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse detonations on the '(+)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for inverse detonations
 *     (on the '(+)' branch)
 */
double alpha_max_inverse_detonation_plus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse detonation solutions on the '(+)' branch.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse detonations on the '(+)' branch.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for inverse detonations
 *     (on the '(+)' branch)
 */
double alpha_min_inverse_detonation_plus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse hybrid solutions.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse hybrids.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for inverse hybrids
 */
double alpha_max_inverse_hybrid(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse hybrid solutions.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse hybrids.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for inverse hybrids
 */
double alpha_min_inverse_hybrid(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse deflagration solutions.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse deflagrations.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for inverse
 *     deflagrations
 */
double alpha_max_inverse_deflagration(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse deflagration solutions.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse deflagrations.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for inverse
 *     deflagrations
 */
double alpha_min_inverse_deflagration(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     mixed direct-inverse hybrid solutions with \f$v_-=c_-\f$.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for mixed hybrids with \f$v_-=c_-\f$.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for mixed hybrids with
 *     \f$v_- = c_-\f$
 */
double alpha_max_mixed_hybrid_minus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     mixed direct-inverse hybrid solutions with \f$v_-=c_-\f$.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for mixed hybrids with \f$v_-=c_-\f$.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for mixed hybrids with
 *     \f$v_- = c_-\f$
 */
double alpha_min_mixed_hybrid_minus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     mixed direct-inverse hybrid solutions with \f$v_+=c_+\f$.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for mixed hybrids with \f$v_+=c_+\f$.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for mixed hybrids with
 *     \f$v_+ = c_+\f$
 */
double alpha_max_mixed_hybrid_plus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     mixed direct-inverse hybrid solutions with \f$v_+=c_+\f$.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for mixed hybrids with \f$v_+=c_+\f$.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for mixed hybrids with
 *     \f$v_+ = c_+\f$
 */
double alpha_min_mixed_hybrid_plus(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for the
 *     given solution type.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for the solution type specified in @p sol_type.
 * 
 * @param sol_type Type of the solution
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength for the specified
 *     solution type
 */
double alpha_max(
    HydroSolutionType sol_type, 
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true, 
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for the
 *     given solution type.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for the solution type specified in @p sol_type.
 * 
 * @param sol_type Type of the solution
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength for the specified
 *     solution type
 */
double alpha_min(
    HydroSolutionType sol_type, 
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true, 
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$.
 * 
 * This function returns the maximal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) over all solution types for the given wall velocity.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{max}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{max}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The maximal value of the transition strength
 */
double alpha_max(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true, 
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$.
 * 
 * This function returns the minimal value of the phase transition strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) over all solution types for the given wall velocity.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_{N,\mathrm{min}}\f$ (\c true, default) or
 *     \f$\alpha_{+,\mathrm{min}}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return The minimal value of the transition strength
 */
double alpha_min(
    double xi_w, 
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    bool return_alpha_N=true, 
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Critical value of \f$c_-\f$ for onset of non-injectivity on the '(+)'
 *     branch of inverse detonations.
 * 
 * Returns the fitted critical sound speed
 * \f[
 *     c_{-,\mathrm{crit}}(c_+) = c_+\, \Bigl[
 *         1 - \bigl(1-c_+^2\bigr)\, P(c_+^2)
 *     \Bigr],
 * \f]
 * where \f$P(x)\f$ is a degree-20 polynomial fit.
 * 
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall
 * @return Critical value \f$c_{-,\mathrm{crit}}(c_+)\f$ at which
 *     non-injectivity sets in
 */
double inverse_detonation_plus_fold_critical_cm(double cp);

/**
 * @brief Check whether \f$(c_+, c_-)\f$ lies in the non-injective sound-speed
 *     region.
 * 
 * Returns \c true if \f$c_+ > c_-\f$ and \f$c_- < c_{-,\mathrm{crit}}(c_+)\f$,
 * where \f$c_{-,\mathrm{crit}}(c_+)\f$ is the fitted boundary from
 * \c profile_calculator::inverse_detonation_plus_fold_critical_cm.
 * 
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall
 * @param cm Speed of sound \f$c_-\f$ behind the wall
 * @return \c true if the sound speeds lie in the region where more than one
 *     solution can exist for some \f$(\xi_w, \alpha_N)\f$ values
 */
bool check_noninjectivity(double cp, double cm);

/**
 * @brief Classify the type of the solution based on \f$\xi_w\f$ and
 *     \f$\alpha_N\f$ in the fluid frame.
 * 
 * Unless \c settings::noninjectivity_check is set to
 * \c settings::CheckMode::SILENT, the sound-speed pair is checked with
 * \c profile_calculator::check_noninjectivity(cp, cm) before branch selection
 * and a warning or error is thrown depending on the value of
 * \c settings::noninjectivity_check.
 * 
 * @param xi_wall Bubble wall velocity \f$\xi_w\f$
 * @param alpha Pseudotrace phase transition strength
 *     \f$\alpha_N\f$ (far ahead of the wall)
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c profile_calculator::profile_solve_mode)
 * @return Type of the hydrodynamic solution
 */
HydroSolutionType classify_wave_fluid_frame(
    double xi_wall, 
    double alpha,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Classify the type of the solution based on \f$v_+\f$ and \f$v_-\f$ in
 *     the wall rest frame and the wall velocity \f$\xi_w\f$.
 * @param vp Fluid velocity \f$v_+\f$ ahead of the wall
 * @param vm Fluid velocity \f$v_-\f$ behind the wall
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Type of the hydrodynamic solution 
 */
HydroSolutionType classify_wave_wall_frame(
    double vp,
    double vm,
    double xi_w,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag
);

/**
 * @brief Compute the fluid velocities \f$v_+\f$ and \f$v_-\f$ in the wall
 *     rest frame.
 * 
 * @param xi_wall Bubble wall velocity \f$\xi_w\f$
 * @param alpha Pseudotrace phase transition strength 
 *     \f$\alpha_+\f$ ahead of the wall
 * @param sol_type Type of the solution
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return The fluid velocities ahead of and behind the wall in the wall rest
 *     frame, \f$v_+\f$ and \f$v_-\f$
 */
std::pair<double, double> compute_fluid_velocities(
    double xi_wall, 
    double alpha, 
    HydroSolutionType sol_type,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag
);

/**
 * @brief Derivative \f$\frac{d\xi}{dv}\f$.
 * 
 * We solve this differential equation to obtain \f$\xi(v)\f$ and then invert to
 * obtain \f$v(\xi)\f$.
 * This is Eq. (2.10) of @ref ref_Barni_2027 "Barni and Madge (2027)".
 * 
 * @param v Fluid velocity \f$v\f$
 * @param xi Self-similar coordinate \f$\xi=r/t\f$
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return Derivative \f$\frac{d\xi}{dv}\f$
 */
double dxi_dv(double v, double xi, double cs=constants::cs_bag);

/**
 * @brief Derivative \f$\frac{dv}{d\xi}\f$.
 * 
 * We solve this differential equation to directly obtain \f$v(\xi)\f$.
 * This is Eq. (2.10) of @ref ref_Barni_2027 "Barni and Madge (2027)".
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$
 * @param v Fluid velocity \f$v\f$
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return Derivative \f$\frac{dv}{d\xi}\f$
 */
double dv_dxi(double xi, double v, double cs=constants::cs_bag);

/**
 * @brief Solve the differential equation \f$\frac{d\xi}{d v}\f$ for the
 *     velocity profile.
 * 
 * This function solves the differential equation using a fourth-order
 * Runge-Kutta method. The first two arguments are references to the vectors in
 * which the solutions are stored. These vectors must already contain the
 * initial condition.
 * 
 * @param xi_vec Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v_vec Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_stop Final fluid velocity \f$v_f\f$
 * @param dv Fluid velocity step size
 * @param stop_at_shock Whether to stop the evolution when a shock front is
 *     encountered (default: \c true)
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return A status code indicating the reason for termination:
 *     0  - reached v_stop
 *     1  - encountered shock 
 *     -1 - entered unphysical range for \f$\xi\f$
 *     -2 - reached super-luminal fluid velocity \f$v\f$
 */
int integrate_xi_of_v(
    std::vector<double>& xi_vec, 
    std::vector<double>& v_vec, 
    double v_stop, 
    double dv, 
    bool stop_at_shock=true,
    double cs=constants::cs_bag
);

/**
 * @brief Solve the differential equation \f$\frac{dv}{d \xi}\f$ for the
 *     velocity profile.
 * 
 * This function solves the differential equation using a fourth-order
 * Runge-Kutta method. The first two arguments are references to the vectors in
 * which the solutions are stored. These vectors must already contain the
 * initial condition.
 * 
 * @param xi_vec Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v_vec Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_stop Final fluid velocity \f$v_f\f$
 * @param dxi Self-similar coordinate step size
 * @param stop_at_shock Whether to stop the evolution when a shock front is
 *     encountered (default: \c true)
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return A status code indicating the reason for termination:
 *     0  - reached v_stop
 *     1  - encountered shock 
 *     -1 - entered unphysical range for \f$\xi\f$
 *     -2 - reached super-luminal fluid velocity \f$v\f$
 *     -3 - crossed Jouguet front
 */
int integrate_v_of_xi(
    std::vector<double>& xi_vec,
    std::vector<double>& v_vec,
    double v_stop,
    double dxi,
    bool stop_at_shock=true,
    double cs=constants::cs_bag
);

/**
 * @brief Solve the differential equation for the velocity profile.
 * 
 * This function solves the differential equation using a fourth-order
 * Runge-Kutta method. The first two arguments are references to the vectors in
 * which the solutions are stored. These vectors must already contain the
 * initial condition.
 * 
 * @param xi_vec Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v_vec Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_stop Final fluid velocity \f$v_f\f$
 * @param step_size Fluid velocity or self-similar coordinate step size
 *     (depending on @p solve_mode)
 * @param stop_at_shock Whether to stop the evolution when a shock front is
 *     encountered (default: \c true)
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return A status code indicating the reason for termination:
 *     0  - reached v_stop
 *     1  - encountered shock 
 *     -1 - entered unphysical range for \f$\xi\f$
 *     -2 - reached super-luminal fluid velocity \f$v\f$
 *     -3 - crossed Jouguet front
 */
int integrate_velocity_profile(
    std::vector<double>& xi_vec,
    std::vector<double>& v_vec,
    double v_stop,
    double step_size,
    bool stop_at_shock=true,
    double cs=constants::cs_bag,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute a rarefaction wave attached to the bubble wall.
 * 
 * This function integrates the rarefaction-wave part of a velocity profile
 * from the bubble wall, given by the last elements of @p xi and @p v, towards
 * the speed of sound \f$\xi = c_s\f$.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$ (must contain the wall position \f$\xi_w\f$ as last element)
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 *     (must contain the fluid velocity \f$v_w\f$ at the wall as last element)
 * @param c_s Speed of sound \f$c_s\f$ in the region of the rarefaction wave
 *     (default: \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default:
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Status code returned by \c integrate_velocity_profile (0 if the
 *     integration reached the final fluid velocity)
 * @throws std::runtime_error if the initial conditions are not valid for a
 *     rarefaction wave
 */
int compute_rarefaction_wave(
    std::vector<double>& xi,
    std::vector<double>& v,
    double c_s=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute a compression wave attached to the bubble wall.
 * 
 * This function integrates the compression-wave part of a velocity profile
 * from the bubble wall, given by the last elements of @p xi and @p v, until a
 * shock front is encountered.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$ (must contain the wall position \f$\xi_w\f$ as last element)
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 *     (must contain the fluid velocity \f$v_w\f$ at the wall as last element)
 * @param c_s Speed of sound \f$c_s\f$ in the region of the compression wave
 *     (default: \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default:
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Status code returned by \c integrate_velocity_profile (1 if a shock
 *     front was encountered)
 * @throws std::runtime_error if the initial conditions are not valid for a
 *     compression wave
 */
int compute_compression_wave(
    std::vector<double>& xi,
    std::vector<double>& v,
    double c_s=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for a direct detonation.
 * 
 * This function solves the fluid velocity profile for a direct detonation and
 * stores the results of the velocity \f$v\f$ as a function of the self-similar
 * coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_direct_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for a direct deflagration.
 * 
 * This function solves the fluid velocity profile for a direct deflagration
 * and stores the results of the velocity \f$v\f$ as a function of the 
 * self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_direct_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for a direct hybrid ('(+/-)' branch).
 * 
 * This function solves the fluid velocity profile for direct hybrids (on 
 * either branch) and stores the results of the velocity \f$v\f$ as a function
 * of the self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_direct_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double xi_wall, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for an inverse detonation ('(+/-)'
 *     branch).
 * 
 * This function solves the fluid velocity profile for an inverse detonation
 * on either branch and stores the results of the velocity \f$v\f$ as a function
 * of the self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_inverse_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for an inverse deflagration.
 * 
 * This function solves the fluid velocity profile for an inverse deflagration
 * and stores the results of the velocity \f$v\f$ as a function of the
 * self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall 
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_inverse_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for an inverse hybrid.
 * 
 * This function solves the fluid velocity profile for an inverse hybrid and
 * stores the results of the velocity \f$v\f$ as a function of the self-similar
 * coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_inverse_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double v_minus, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for a mixed hybrid with \f$v_-=c_-\f$.
 * 
 * This function solves the fluid velocity profile for mixed hybrids with
 * \f$v_-=c_-\f$ and stores the results of the velocity \f$v\f$ as a function
 * of the self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_mixed_hybrid_minus(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double xi_wall, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for a mixed hybrid with \f$v_+=c_+\f$.
 * 
 * This function solves the fluid velocity profile for mixed hybrids with
 * \f$v_+=c_+\f$ and stores the results of the velocity \f$v\f$ as a function
 * of the self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param c_plus Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param c_minus Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_mixed_hybrid_plus(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double v_minus, 
    double c_plus = constants::cs_bag,
    double c_minus = constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Get the position of the shock front.
 * 
 * @param sol_type Type of the solution
 * @param size Length of the profile vectors
 * @return Index of the first element outside the shock front
 */
size_t get_shock_index(HydroSolutionType sol_type, size_t size);

/**
 * @brief Get the coordinate and velocity of the shock front.
 * 
 * @param sol_type Type of the solution
 * @param xi Self-similar coordinate \f$\xi=r/t\f$ of the profile
 * @param v Fluid velocity profile \f$v(\xi)\f$
 * @return The shock front position \f$(\xi_\text{sh}, v_\text{sh})\f$
 */
std::pair<double, double> get_shock_front(
    HydroSolutionType sol_type, 
    const std::vector<double>& xi, 
    const std::vector<double>& v
);

/**
 * @brief Compute the velocity profile for a bubble in a cosmological phase
 *     transition.
 * 
 * This function solves the fluid velocity profile for a cosmological phase
 * transition and stores the results of the velocity \f$v\f$ as a function of
 * the self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param alpha_N Pseudotrace transition strength \f$\alpha_N\f$ (far in front
 *     of the wall)
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall 
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_N, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the velocity profile for a bubble in a cosmological phase
 *     transition.
 * 
 * This function solves the fluid velocity profile for a cosmological phase
 * transition and stores the results of the velocity \f$v\f$ as a function of
 * the self-similar coordinate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * If a particular \f$\alpha_+\f$ root is already known, this overload
 * constructs the corresponding stationary profile directly for the selected
 * @p sol_type.
 * 
 * @param sol_type Type of the solution
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param alpha_p Pseudotrace transition strength \f$\alpha_+\f$ ahead of the
 *     wall
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Index of the first element outside the bubble (i.e. at the wall 
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile(
    HydroSolutionType sol_type, 
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_p, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Calculate a segment of the enthalpy density integral.
 * 
 * This function calculates the segment of the integral in Eq. (2.12) of
 * @ref ref_Barni_2027 "Barni and Madge (2027)"
 * \f[
 *  I = \left(1+\frac{1}{c_s^2}\right) \int\limits_{v_1}^{v_2} \gamma^2 \mu\, dv
 *    = \left(1+\frac{1}{c_s^2}\right) \int\limits_{v_1}^{v_2} \frac{1}{1-v^2}
 *      \frac{\xi - v}{1 - \xi v}\, dv
 * \f]
 * assuming that \f$\xi\f$ is linear in \f$v\f$ within the segment
 * \f$(\xi_1,\xi_2)\f$, i.e.
 * \f$\xi = \xi_0 + \partial_v \xi\, v\f$, where
 * \f$\xi_0=\frac{\xi_1 v_2 - \xi_2 v_1}{v_2-v_1}\f$ and
 * \f$\partial_v \xi = \frac{\xi_2-\xi_1}{v_2-v_1}\f$.
 * The integral can then be calculated using partial fraction decomposition:
 * \f[
 *  I = -\left(1+\frac{1}{c_s^2}\right) \int\limits_{v_1}^{v_2} \left[
 *     \frac{v}{1-v^2}
 *     + \frac{\partial_v \xi\,v + \xi_0}{\partial_v \xi\, v^2 + \xi_0 v - 1}
 *   \right] dv
 *    = \left(1+\frac{1}{c_s^2}\right) \left\{
 *          \frac{1}{2} \int\limits_{1-v_1^2}^{1-v_2^2} \frac{du}{u}
 *          + \frac{1}{2} \int\limits_{1-\xi_2 v_2}^{1-\xi_1 v_1} \frac{du}{u}
 *          - \frac{1}{2} \int\limits_{v_1}^{v_2}
 *              \frac{\xi_0}{\partial_v \xi\, v^2 + \xi_0 v - 1} dv
 *       \right\} .
 * \f]
 * 
 * @param xi1 Self-similar coordinate \f$\xi_1\f$ at the beginning of the 
 *     segment
 * @param xi2 Self-similar coordinate \f$\xi_2\f$ at the end of the segment
 * @param v1 Fluid velocity \f$v_1\f$ at the beginning of the segment
 * @param v2 Fluid velocity \f$v_2\f$ at the end of the segment
 * @param cs Speed of sound \f$c_s\f$ (default: \c constants::cs_bag)
 * @return Value of the integral \f$I\f$ for the segment
 */
double integrate_enthalpy_segment(
    double xi1,
    double xi2,
    double v1,
    double v2,
    double cs = constants::cs_bag
);

/**
 * @brief Compute the enthalpy density profile.
 * 
 * Compute the enthalpy density profile from the fluid velocity profile,
 * Eq. (2.12) of @ref ref_Barni_2027 "Barni and Madge (2027)".
 * 
 * @param xi Reference to the vector containing the self-similar coordinate 
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param idx_shock Index of the first element outside the shock front (use
 *     \c constants::NO_INDEX if there is no shock front, default:
 *     \c constants::NO_INDEX)
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Enthalpy density profile \f$w(\xi)\f$
 */
std::vector<double> compute_enthalpy_density(
    const std::vector<double>& xi, 
    const std::vector<double>& v, 
    size_t idx_wall, 
    size_t idx_shock=constants::NO_INDEX,
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the energy density profile.
 * 
 * We assume the constant-sound-speed (or \f$\mu\f$-\f$\nu\f$) model,
 * cf. @ref ref_Leitao_2015 "Leitao and Megevand (2015)" and
 * @ref ref_Giese_2020 "Giese et al. (2020)", i.e.
 * \f$e(\xi) = \frac{1}{1+c_s^2(\xi)} w(\xi) + \epsilon(\xi)\f$, where
 * \f$c_s(\xi) = c_+\Theta(\xi-\xi_w) + c_-\Theta(\xi_w-\xi)\f$ and
 * \f$\epsilon(\xi)=\epsilon_+\Theta(\xi-\xi_w)+\epsilon_-\Theta(\xi_w-\xi)\f$
 * is the vacuum energy contribution. Given the bag constant
 * \f$\epsilon = \epsilon_+ - \epsilon_-\f$ (@p epsilon), we set
 * \f$\epsilon_+=\epsilon\f$, \f$\epsilon_-=0\f$ if \f$\epsilon \geq 0\f$, and
 * \f$\epsilon_+ = 0\f$, \f$\epsilon_- = -\epsilon\f$ otherwise.
 * 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param epsilon Bag constant \f$\epsilon = \epsilon_+ - \epsilon_-\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Energy density profile \f$e(\xi)\f$
 */
std::vector<double> compute_energy_density_profile(
    const std::vector<double>& w,
    size_t idx_wall,
    double epsilon,
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the pressure profile.
 * 
 * We assume the constant-sound-speed (or \f$\mu\f$-\f$\nu\f$) model,
 * cf. @ref ref_Leitao_2015 "Leitao and Megevand (2015)" and
 * @ref ref_Giese_2020 "Giese et al. (2020)", i.e.
 * \f$p(\xi) = \frac{c_s^2(\xi)}{1+c_s^2(\xi)} w(\xi) - \epsilon(\xi)\f$, where
 * \f$c_s(\xi) = c_+\Theta(\xi-\xi_w) + c_-\Theta(\xi_w-\xi)\f$ and
 * \f$\epsilon(\xi)=\epsilon_+\Theta(\xi-\xi_w)+\epsilon_-\Theta(\xi_w-\xi)\f$
 * is the vacuum energy contribution. Given the bag constant
 * \f$\epsilon = \epsilon_+ - \epsilon_-\f$ (@p epsilon), we set
 * \f$\epsilon_+=\epsilon\f$, \f$\epsilon_-=0\f$ if \f$\epsilon \geq 0\f$, and
 * \f$\epsilon_+ = 0\f$, \f$\epsilon_- = -\epsilon\f$ otherwise.
 * 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param epsilon Bag constant \f$\epsilon = \epsilon_+ - \epsilon_-\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Pressure profile \f$p(\xi)\f$
 */
std::vector<double> compute_pressure_profile(
    const std::vector<double>& w, 
    size_t idx_wall, 
    double epsilon,
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the temperature profile.
 * 
 * We assume the constant-sound-speed (or \f$\mu\f$-\f$\nu\f$) model,
 * cf. @ref ref_Leitao_2015 "Leitao and Megevand (2015)" and
 * @ref ref_Giese_2020 "Giese et al. (2020)", i.e. \f$w = (1+c_s^2) a T^\nu\f$,
 * where \f$\nu = 1 + \frac{1}{c_s^2}\f$ and \f$a\f$ is proportional to the
 * effective number of relativistic degrees of freedom in the phase. We use
 * the equal-temperature enthalpy ratio across the wall at the nucleation
 * temperature, \f$\Psi_N = \frac{w_-(T_N)}{w_+(T_N)}\f$, to match the
 * temperatures at the bubble wall.
 * 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param T_N Nucleation temperature \f$T_N\f$ far ahead of the wall
 * @param Psi_N Enthalpy ratio \f$\Psi_N\f$ across the wall at nucleation
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Temperature profile \f$T(\xi)\f$
 */
std::vector<double> compute_temperature_profile(
    const std::vector<double>& w,
    size_t idx_wall,
    double T_N,
    double Psi_N,
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the entropy density profile.
 * 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param T Reference to the vector containing the temperature \f$T\f$
 * @return Entropy density profile \f$s(\xi)\f$
 */
std::vector<double> compute_entropy_density_profile(
    const std::vector<double>& w,
    const std::vector<double>& T
);

/**
 * @brief Check the entropy flux at discontinuities of a profile.
 * 
 * Check that the entropy flux at a discontinuity does not violate the second
 * law of thermodynamics, i.e. that it does not decrease in the rest frame of
 * the discontinuity, \f$s_- \gamma_- v_- \geq s_+ \gamma_+ v_+\f$.
 * 
 * @note The entropy density depends on the choice of \f$\Psi_N\f$, cf.
 *     \c compute_temperature_profile.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param s Reference to the vector containing the entropy density \f$s\f$
 * @return Index of the first element outside the first discontinuity that
 *     violates the condition, or \c constants::NO_INDEX if the profile is fine
 */
size_t check_entropy_flux(
    const std::vector<double>& xi,
    const std::vector<double>& v,
    const std::vector<double>& s
);

/**
 * @brief Compute the bag strength parameter \f$\alpha_{\mathrm{bag}}\f$ from
 *     the pseudotrace strength parameter \f$\alpha_{\bar\theta}\f$.
 * 
 * @param alpha Pseudotrace transition strength \f$\alpha_{\bar\theta}\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Bag strength parameter \f$\alpha_{\mathrm{bag}}\f$
 */
double compute_alpha_bag_from_pseudotrace(
    double alpha, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the pseudotrace strength parameter \f$\alpha_{\bar\theta}\f$
 *     from the bag strength parameter \f$\alpha_{\mathrm{bag}}\f$.
 * 
 * @param alpha Bag strength parameter \f$\alpha_{\mathrm{bag}}\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Pseudotrace transition strength \f$\alpha_{\bar\theta}\f$
 */
double compute_alpha_pseudotrace_from_bag(
    double alpha, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the pseudotrace strength parameter at some point (in front of
 *     the wall) from the value at some other point and the enthalpy ratio.
 * 
 * @param alpha Pseudotrace transition strength \f$\alpha_i\f$
 * @param r_w Enthalpy ratio \f$r_w = \frac{w_i}{w_f}\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Pseudotrace transition strength \f$\alpha_f\f$
 */
double rescale_alpha_by_enthalpy(
    double alpha, 
    double r_w,
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the bag constant (vacuum energy density).
 * 
 * \f[
 *     \epsilon = \epsilon_+ - \epsilon_- = \frac{3\,c_-^2\, w_+}{1+c_-^2}
 *     \left(
 *          \alpha_+ - \frac{c_-^2-c_+^2}{3\,c_-^2 (1+c_+^2)}
 *     \right) = \frac{3\,c_-^2\, w_N}{1+c_-^2} \left(
 *          \alpha_N - \frac{c_-^2-c_+^2}{3\,c_-^2 (1+c_+^2)}
 *     \right)
 * \f]
 * 
 * @param w Enthalpy density, either \f$w_+\f$ (directly in front) or \f$w_N\f$
 *     (far in front of the bubble wall)
 * @param alpha Pseudotrace transition strength, either \f$\alpha_+\f$ (directly
 *     in front) or \f$\alpha_N\f$ (far in front of the bubble wall)
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @return Bag constant \f$\epsilon = \epsilon_+ - \epsilon_-\f$
 */
double compute_epsilon(
    double w, 
    double alpha, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag
);

/**
 * @brief Compute the pseudotrace transition strength parameter \f$\alpha_N\f$
 *     from the value in front of the bubble wall \f$\alpha_+\f$.
 * 
 * @param sol_type Type of the solution
 * @param alpha_plus Pseudotrace transition strength parameter \f$\alpha_+\f$
 *     immediately in front of the bubble wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Pseudotrace transition strength \f$\alpha_N\f$ far in front of the
 *     bubble wall
 */
double compute_alpha_N(
    HydroSolutionType sol_type, 
    double alpha_plus, 
    double xi_wall, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Find the pseudotrace transition strength \f$\alpha_+\f$ from the value
 *     \f$\alpha_N\f$ in the false vacuum phase.
 * 
 * The solution type is determined by @ref classify_wave_fluid_frame. The 
 * automatic classification may report a non-injective sound-speed region
 * depending on the \c settings::noninjectivity_check switch.
 * 
 * The implementation performs a single-root inversion and returns one
 * \f$\alpha_+\f$. On a folded branch, more than one \f$\alpha_+\f$ can
 * correspond to the same \f$\alpha_N\f$ even within the same hydrodynamic
 * branch.
 * 
 * @param alpha_N Pseudotrace transition strength \f$\alpha_N\f$ far in front of
 *     the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return One value of \f$\alpha_+\f$
 */
double find_alpha_plus(
    double alpha_N, 
    double xi_wall, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Find the pseudotrace transition strength \f$\alpha_+\f$ from the value
 *     \f$\alpha_N\f$ in the false vacuum phase.
 * 
 * The search is restricted to the selected @p sol_type. The implementation
 * performs a single-root inversion and returns one \f$\alpha_+\f$. On a folded
 * branch, more than one \f$\alpha_+\f$ can correspond to the same
 * \f$\alpha_N\f$ even within the same @p sol_type. This function will, however,
 * find at most one root.
 * 
 * @param sol_type Type of the solution
 * @param alpha_N Pseudotrace transition strength \f$\alpha_N\f$ far in front of
 *     the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return One value of \f$\alpha_+\f$
 */
double find_alpha_plus(
    HydroSolutionType sol_type,
    double alpha_N, 
    double xi_wall, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Compute the efficiency.
 * 
 * Compute the kinetic energy fraction
 * \f$K = \frac{\rho_\mathrm{kin}}{\rho_\mathrm{tot}}\f$
 * or the efficiency factor \f$\kappa_\mathrm{sw}\f$ for converting vacuum
 * energy to bulk kinetic energy, Eq. (6.3) and (6.4) of
 * @ref ref_Barni_2027 "Barni and Madge (2027)". See also Eq. (2.30) of
 * @ref ref_Espinosa_2010 "Espinosa and No (2010)" or Eq. (4.10) of
 * @ref ref_Barni_2024 "Barni et al. (2024)".
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$ 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param alpha_N Pseudotrace transition strength \f$\alpha_N\f$ (far in front
 *     of the wall)
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param normalization \c EfficiencyNormalization switch to determine how to
 *     normalize the efficiency (default:
 *     \c EfficiencyNormalization::TotalEnergyDensity)
 * @return Efficiency \f$\kappa_\mathrm{sw}\f$ or \f$K\f$ (depending on
 *     @p normalization)
 */
double compute_efficiency(
    const std::vector<double>& xi,
    const std::vector<double>& v,
    const std::vector<double>& w,
    size_t idx_wall,
    double alpha_N,
    double cp = constants::cs_bag,
    double cm = constants::cs_bag,
    EfficiencyNormalization normalization 
        = EfficiencyNormalization::TotalEnergyDensity
);

/**
 * @brief Compute the enthalpy ratio \f$\Psi_+\f$ in local thermal equilibrium.
 * 
 * Compute the enthalpy ratio \f$\Psi_+ = \frac{w_-(T_+)}{w_+(T_+)}\f$ ahead
 * of and behind the wall, evaluated at the temperature ahead of the wall.
 * This assumes local thermal equilibrium and uses entropy conservation to
 * relate the temperatures on both sides of the wall, i.e.
 * \f$s_+ v_+ \gamma_+ = s_- v_- \gamma_-\f$.
 * 
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param alpha_N Pseudotrace transition strength \f$\alpha_N\f$ (far in front
 *     of the wall)
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Enthalpy ratio \f$\Psi_+\f$ across the wall
 */
double compute_enthalpy_ratio_LTE(
    double xi_wall, 
    double alpha_N, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/**
 * @brief Find the LTE wall velocity \f$\xi_w\f$ for a given enthalpy ratio
 *     \f$\Psi_+\f$.
 * 
 * Find the wall velocity \f$\xi_w\f$ for a given enthalpy ratio
 * \f$\Psi_+ = \frac{w_-(T_+)}{w_+(T_+)}\f$ assuming local thermal equilibrium.
 * 
 * @param alpha_N Pseudotrace transition strength \f$\alpha_N\f$ (far in front
 *     of the wall)
 * @param Psi Enthalpy ratio \f$\Psi_+\f$
 * @param xi_min Minimal wall velocity \f$\xi_w^\mathrm{min}\f$ to consider
 * @param xi_max Maximal wall velocity \f$\xi_w^\mathrm{max}\f$ to consider
 * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
 *     \c constants::cs_bag)
 * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
 *     \c constants::cs_bag)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 
 *     \c profile_calculator::profile_step_size)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c profile_calculator::profile_solve_mode)
 * @return Wall velocity \f$\xi_w\f$
 */
double find_LTE_wall_velocity(
    double alpha_N, 
    double Psi, 
    double xi_min, 
    double xi_max, 
    double cp=constants::cs_bag,
    double cm=constants::cs_bag,
    double step_size=profile_step_size,
    ProfileSolveMode solve_mode=profile_solve_mode
);

/** 
 * @brief Check profile for consistency.
 * 
 * This function returns @c true if @p xi and @p profile have the same size and
 * at least size 6. Otherwise it throws an error (if @p throw_error is @c true)
 * or returns @c false.
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$ of the profile
 * @param profile Vector containing the fluid profile
 * @param throw_error Whether to throw an error if the test fails (default:
 *     \c true)
 * @return Whether the profile passes the test
 */
bool check_profiles(
    const std::vector<double>& xi, 
    const std::vector<double>& profile, 
    bool throw_error=true
);


} // namespace profile_calculator