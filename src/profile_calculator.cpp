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

#include "hydrodex/profile_calculator.hpp"

#include "hydrodex/constants.hpp"
#include "hydrodex/settings.hpp"
#include "hydrodex/utils.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/math/tools/roots.hpp>

using namespace hydrodex::constants;
using namespace hydrodex::settings;

namespace hydrodex::profile_calculator {

double equal_rel_tol = 1e-9;
double equal_abs_tol = 1e-12;
double profile_step_size = 1e-4;
ProfileSolveMode profile_solve_mode = ProfileSolveMode::D_XI_D_V;
double alpha_rel_tol = 1e-6; 
double alpha_abs_tol = 1e-9; 
double xi_wall_rel_tol = 1e-6;
double xi_wall_abs_tol = 1e-9;

const std::string to_string(HydroSolutionType type) {
    switch (type) {
        case HydroSolutionType::Deflagration:        
            return "Direct Deflagration";
        case HydroSolutionType::Detonation:          
            return "Direct Detonation";
        case HydroSolutionType::Hybrid:              
            return "Direct Hybrid ('-' branch)";
        case HydroSolutionType::HybridPlus:              
            return "Direct Hybrid ('+' branch)";
        case HydroSolutionType::InverseDeflagration: 
            return "Inverse Deflagration";
        case HydroSolutionType::InverseDetonation:   
            return "Inverse Detonation ('-' branch)";
        case HydroSolutionType::InverseDetonationPlus:   
            return "Inverse Detonation ('+' branch)";
        case HydroSolutionType::InverseHybrid:       
            return "Inverse Hybrid";
        case HydroSolutionType::MixedHybridMinus:
            return "Mixed Hybrid (v- = c-)";
        case HydroSolutionType::MixedHybridPlus:
            return "Mixed Hybrid (v+ = c+)";
        default:                                     
            return "Unknown";
    }
}

std::ostream& operator<<(std::ostream& os, HydroSolutionType type) { 
    return os << to_string(type); 
}

bool are_equal(double x, double y, double eps_rel, double eps_abs) {
    double precision = std::max(
        eps_abs, 
        eps_rel * std::max(std::abs(x), std::abs(y))
    );
    return std::abs(x - y) < precision;
}

double shift_safe(double x, double y, double eps_rel, double eps_abs) {
    double eps = std::min(eps_abs, eps_rel * std::abs(x - y));
    return (x < y ? x + eps : x - eps);
}

double gamma(double v) { return 1.0 / std::sqrt(1.0 - v * v); }

double mu(double xi, double v) { return (xi - v) / (1.0 - v * xi); }

bool check_velocity(double v) { return (v > 0.0 && v < 1.0); }

double get_nu_from_cs(double cs) {
    if ( ! check_velocity(cs) ) { throw std::invalid_argument(
        "Invalid sound speed: cs = " + std::to_string(cs) + 
        ". Must be 0 < cs < 1."
    ); }
    return 1.0 + 1.0 / (cs * cs);
}

double get_cs_from_nu(double nu) {
    if (nu <= 2.0) {
        throw std::invalid_argument(
            "Invalid nu parameter: nu = " + std::to_string(nu) +
            ". Must be nu > 2."
        );
    }
    return 1.0 / std::sqrt(nu - 1.0);
}

double shock_condition(double xi, double v, double cs) { 
    return mu(xi, v) * xi - cs * cs; 
}

double jouguet_condition(double xi, double v, double cs) {
    return mu(xi, v) * mu(xi, v) - cs * cs;
}

std::pair<double, double> get_v_plus_from_v_minus(
    double vm,
    double alpha,
    double cm
) {
    double cm2 = cm * cm;
    double term = 0.5 * (vm + cm2 / vm);
    double disc = std::sqrt(
        term * term
        - cm2 * (1.0 - 3.0 * alpha) * (1.0 + 3.0 * cm2 * alpha)
    );
    double denom = 1.0 + 3.0 * cm2 * alpha;
    double vplus1 = (term - disc) / denom;
    double vplus2 = (term + disc) / denom;
    // handle limits alpha -> - 1/(3 cm^2):
    // '-' branch is finite, '+' branch is infinite
    if ( ! std::isfinite(1.0 / denom) ) { 
        return {
            0.5 * (1.0 + cm2) / term, 
            2.0 * term / denom - 0.5 * (1.0 + cm2) / term
        }; 
    }
    return {vplus1, vplus2};
}

double get_v_plus_from_v_minus(
    double vm, 
    double alpha,
    RootBranch branch,
    double cm
) {
    double cm2 = cm * cm;
    double term = 0.5 * (vm + cm2 / vm);
    double disc = std::sqrt(
        term * term
        - cm2 * (1.0 - 3.0 * alpha) * (1.0 + 3.0 * cm2 * alpha)
    );
    double denom = 1.0 + 3.0 * cm2 * alpha;
    if ( ! std::isfinite(1.0 / denom) ) { 
        if ( branch == RootBranch::Minus ) { 
            return 0.5 * (1.0 + cm2) / term; 
        } else { return 2.0 * term / denom - 0.5 * (1.0 + cm2) / term; }
    }
    return (term + branch * disc) / denom;
}

std::pair<double, double> get_v_minus_from_v_plus(
    double vp,
    double alpha,
    double cm   
) {
    double cm2 = cm * cm;
    double term1 = 0.5 * (1.0 + 3.0 * cm2 * alpha) * vp;
    double term2 = 0.5 * cm2 * (1.0 - 3.0 * alpha) / vp;
    double sum = term1 + term2;
    double disc = std::sqrt(sum * sum - cm2);
    double vminus1 = sum - disc;
    double vminus2 = sum + disc;
    return {vminus1, vminus2};
}

double get_v_minus_from_v_plus(
    double vp,
    double alpha,
    RootBranch branch,
    double cm   
) {
    double cm2 = cm * cm;
    double term1 = 0.5 * (1.0 + 3.0 * cm2 * alpha) * vp;
    double term2 = 0.5 * cm2 * (1.0 - 3.0 * alpha) / vp;
    double sum = term1 + term2;
    double disc = std::sqrt(sum * sum - cm2);
    return sum + branch * disc;
}

double compute_alpha_plus(double vp, double vm, double cm) {
    double factor1 = (vp * vm - cm * cm); 
    double factor2 = (vp - vm);
    double denominator = 3.0 * cm * cm * (1.0 - vp * vp) * vm;
    return factor1 * factor2 / denominator;
}

double get_v_plus_extremum_from_v_minus(double vm, double cm) {
    double cm2 = cm * cm;
    if ( vm <= cm2 ) { return 1.0; }
    double numerator = cm2 + vm * vm;
    double term1 = (1.0 + cm2) * vm;
    double term2 = std::sqrt( (1.0 - vm * vm ) * (vm * vm - cm2 * cm2) );
    return (term1 - term2) / numerator;
}

double get_v_plus_extremum_from_alpha(double alpha, double cm) {
    if ( alpha > 0.0 ) {
        throw std::invalid_argument(
            "Attempted to calculate extremum position of v-(v+, alpha) for"
            " positive alpha: " + std::to_string(alpha) + "."
        );
    }
    if ( alpha < (cm * cm - 1.0)/(6.0 * cm * cm) ) {
        return 1.0;
    }
    double numerator = 1.0 - 3.0 * alpha;
    double denominator = 1.0 + 3.0 * cm * cm * alpha;
    return cm * std::sqrt(numerator / denominator);
}

double get_v_minus_extremum_from_v_plus(
    double vp, 
    RootBranch branch, 
    double cm
) {
    if ( vp < cm ) {
        throw std::invalid_argument(
            "Attempted to calculate extremum of v- for v+ < c-. Got v+ = "
            + std::to_string(vp) + " and c- = " + std::to_string(cm) + "."
        );
    }
    double cm2 = cm * cm;
    double denominator = 1.0 + vp * vp;
    double term1 = (1.0 + cm2) * vp;
    double term2 = std::sqrt( (vp * vp - cm2) * (1.0 - cm2 * vp * vp) );
    return (term1 + branch * term2) / denominator;
}

double get_v_minus_extremum_from_alpha(
    double alpha, 
    RootBranch branch,
    double cm
) {
    if ( alpha > 0.0 ) {
        throw std::invalid_argument(
            "Attempted to calculate extremum position of v-(v+, alpha) for"
            " positive alpha: " + std::to_string(alpha) + "."
        );
    }
    if ( alpha < (cm * cm - 1.0)/(6.0 * cm * cm) ) {
        switch ( branch ) {
            case RootBranch::Minus:
                return cm * cm;
            case RootBranch::Plus:
                return 1.0;
        }
    }
    double factor1 = 1.0 + 3.0 * cm * cm * alpha;
    double factor2 = 1.0 - 3.0 * alpha;
    double sqrt1 = std::sqrt(factor1 * factor2);
    double sqrt2 = std::sqrt(factor1 * factor2 - 1.0);
    return cm * (sqrt1 + branch * sqrt2);
}

double direct_jouguet_velocity(double alpha, double cm) {
    if ( alpha < 0.0 ) { throw std::invalid_argument(
        "Attempted to calculate direct Jouguet velocity for negative alpha: " + 
        std::to_string(alpha) + "."
    ); }
    return get_v_plus_from_v_minus(cm, alpha, RootBranch::Plus, cm);
}

double inverse_jouguet_velocity(double alpha, double cp, double cm) {
    if ( alpha > 0.0 ) { throw std::invalid_argument(
        "Attempted to calculate inverse Jouguet velocity for positive alpha: " + 
        std::to_string(alpha) + "."
    ); }
    return std::min(
        get_v_minus_from_v_plus(cp, alpha, RootBranch::Minus, cm),
        get_v_minus_extremum_from_alpha(alpha, RootBranch::Minus, cm)
    );
}

double xi_min(HydroSolutionType sol_type, double cp, double cm) {
    switch(sol_type) {
        case HydroSolutionType::Detonation:
            return std::max(cp, cm);
        case HydroSolutionType::HybridPlus:
        case HydroSolutionType::MixedHybridMinus:
        case HydroSolutionType::MixedHybridPlus:
            return ( cm < cp ? cm : 1.0 );
        case HydroSolutionType::Hybrid:
            return cm;
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::InverseDetonation:
            return 0.0;
        case HydroSolutionType::InverseDetonationPlus:
            return ( cm < cp ? get_v_minus_extremum_from_v_plus(
                cp, RootBranch::Minus, cm
            ): 1.0 );
        case HydroSolutionType::InverseHybrid:
            return ( cm * cm < cp ? cm * cm : 1.0 );
        case HydroSolutionType::InverseDeflagration:
            return std::max(cm * cm, cp);
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
}

double xi_max(HydroSolutionType sol_type, double cp, double cm) {
    switch(sol_type) {
        case HydroSolutionType::Detonation:
        case HydroSolutionType::Hybrid:
        case HydroSolutionType::InverseDeflagration:
            return 1.0;
        case HydroSolutionType::HybridPlus:
            return ( cm < cp ? std::min(cp * cp / cm, 1.0) : 0.0 );
        case HydroSolutionType::Deflagration:
            return cm;
        case HydroSolutionType::InverseDetonation:
            return std::min(cp, cm);
        case HydroSolutionType::InverseDetonationPlus:
            return ( cm < cp ? cm : 0.0 );
        case HydroSolutionType::InverseHybrid:
            return ( cm * cm < cp ? cp : 0.0 );
        case HydroSolutionType::MixedHybridMinus:
        case HydroSolutionType::MixedHybridPlus:
            return ( cm < cp ? cp : 0.0 );
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
}

void check_wall_velocity(
    HydroSolutionType sol_type,
    double xi_w,
    double cp,
    double cm
) {
    double ximin = xi_min(sol_type, cp, cm);
    double ximax = xi_max(sol_type, cp, cm);
    if ( ximin > ximax ) { throw std::invalid_argument(
        "Solutions of type " + to_string(sol_type) + " do not exist for c+ = "
        + std::to_string(cp) + " and c- = " + std::to_string(cm) +"."
    ); }
    if ( xi_w < ximin || xi_w > ximax ) { throw std::invalid_argument(
        "For c+ = " + std::to_string(cp) + " and c- = " + std::to_string(cm)
        + ", the wall velocity for solutions of type " + to_string(sol_type)
        + " must be within " + std::to_string(ximin) + " < xi < "
        + std::to_string(ximax) + ". Got xi = " + std::to_string(xi_w) + "."
    ); }
}

double alpha_max_detonation(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::Detonation, xi_w, cp, cm);
    return compute_alpha_plus(xi_w, cm , cm);
}

double alpha_min_detonation(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::Detonation, xi_w, cp, cm);
    return 0.0;
}

double alpha_max_hybrid_plus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::HybridPlus, xi_w, cp, cm);
    double vp = std::min(xi_w, cp * cp / xi_w);
    double alpha_p = compute_alpha_plus(vp, cm, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::HybridPlus, 
            shift_safe(0.0, alpha_p, alpha_rel_tol, alpha_abs_tol), 
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_min_hybrid_plus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::HybridPlus, xi_w, cp, cm);
    if ( return_alpha_N ) {
        return compute_alpha_plus(xi_w, cm, cm);
    }
    return 0.0;
}

double alpha_max_hybrid(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::Hybrid, xi_w, cp, cm);
    double alpha_p = 1.0 / 3.0;
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::Hybrid, 
            shift_safe(alpha_p, 0.0, alpha_rel_tol, alpha_abs_tol), 
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_min_hybrid(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::Hybrid, xi_w, cp, cm);
    if ( xi_w < cp * cp / cm ) {
        if ( return_alpha_N ) {
            return compute_alpha_N(
                HydroSolutionType::Hybrid, 
                shift_safe(0.0, 1.0 / 3.0, alpha_rel_tol, alpha_abs_tol), 
                xi_w, 
                cp,
                cm,
                step_size,
                solve_mode
            );
        } 
        return 0.0; 
    }
    if ( return_alpha_N ) {
        return compute_alpha_plus(xi_w, cm , cm);
    }
    return compute_alpha_plus(cp * cp / xi_w, cm, cm);
}

double alpha_max_deflagration(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::Deflagration, xi_w, cp, cm);
    double alpha_p = 1.0 / 3.0;
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::Deflagration, 
            shift_safe(
                alpha_p, 
                alpha_min_deflagration(xi_w, cp, cm, false),
                alpha_rel_tol, 
                alpha_abs_tol
            ), 
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_min_deflagration(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::Deflagration, xi_w, cp, cm);
    if ( return_alpha_N ) { return 0.0; }
    if ( xi_w > cp ) { return compute_alpha_plus(cp * cp / xi_w, xi_w, cm); }
    return 0.0;
}

double alpha_max_inverse_detonation(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseDetonation, xi_w, cp, cm);
    return 0.0;
}

double alpha_min_inverse_detonation(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseDetonation, xi_w, cp, cm);
    double vp = cp;
    if ( cp > cm ) {
        vp = std::min(vp, get_v_plus_extremum_from_v_minus(xi_w, cm));
    }
    double alpha_p = compute_alpha_plus(vp, xi_w, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::InverseDetonation, 
            shift_safe(alpha_p, 0.0, alpha_rel_tol, alpha_abs_tol), 
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_max_inverse_detonation_plus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseDetonationPlus, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(cp, xi_w, cm);
    double vp_extr = get_v_plus_extremum_from_v_minus(xi_w, cm);
    if ( return_alpha_N ){
        return compute_alpha_N(
            HydroSolutionType::InverseDetonationPlus, 
            shift_safe(
                compute_alpha_plus(vp_extr, xi_w, cm), 
                alpha_p,
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_min_inverse_detonation_plus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseDetonationPlus, xi_w, cp, cm);
    double vp_extr = get_v_plus_extremum_from_v_minus(xi_w, cm);
    double alpha_p = compute_alpha_plus(vp_extr, xi_w, cm);
    if ( return_alpha_N ){
        return compute_alpha_N(
            HydroSolutionType::InverseDetonationPlus, 
            shift_safe(
                compute_alpha_plus(cp, xi_w, cm), 
                alpha_p,
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_max_inverse_hybrid(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseHybrid, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(cp, xi_w, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::InverseHybrid, 
            shift_safe(
                alpha_p, 
                compute_alpha_plus(cp, 1.0, cm),
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_min_inverse_hybrid(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseHybrid, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(cp, 1.0, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::InverseHybrid, 
            shift_safe(
                alpha_p, 
                compute_alpha_plus(cp, xi_w, cm),
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p;
}

double alpha_max_inverse_deflagration(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseDeflagration, xi_w, cp, cm);
    return 0.0;
}

double alpha_min_inverse_deflagration(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::InverseDeflagration, xi_w, cp, cm);
    return compute_alpha_plus(xi_w, 1.0, cm);
}

double alpha_max_mixed_hybrid_minus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::MixedHybridMinus, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(cp, cm, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::MixedHybridMinus, 
            shift_safe(
                compute_alpha_plus(xi_w, cm, cm),
                alpha_p, 
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p; 
}

double alpha_min_mixed_hybrid_minus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::MixedHybridMinus, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(xi_w, cm, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::MixedHybridMinus, 
            shift_safe(
                compute_alpha_plus(cp, cm, cm),
                alpha_p, 
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p; 
}

double alpha_max_mixed_hybrid_plus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::MixedHybridPlus, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(cp, cm, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::MixedHybridPlus, 
            shift_safe(
                alpha_p, 
                compute_alpha_plus(cp, xi_w, cm),
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p; 
}

double alpha_min_mixed_hybrid_plus(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    check_wall_velocity(HydroSolutionType::MixedHybridPlus, xi_w, cp, cm);
    double alpha_p = compute_alpha_plus(cp, xi_w, cm);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::MixedHybridPlus, 
            shift_safe(
                alpha_p, 
                compute_alpha_plus(cp, cm, cm),
                alpha_rel_tol, 
                alpha_abs_tol
            ),
            xi_w, 
            cp,
            cm,
            step_size,
            solve_mode
        );
    }
    return alpha_p; 
}

double alpha_max(
    HydroSolutionType sol_type, 
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    switch ( sol_type ) {
        case HydroSolutionType::Detonation:
            return alpha_max_detonation(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
       case HydroSolutionType::HybridPlus: 
            return alpha_max_hybrid_plus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
       case HydroSolutionType::Hybrid: 
            return alpha_max_hybrid(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::Deflagration: 
            return alpha_max_deflagration(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDetonation:
            return alpha_max_inverse_detonation(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDetonationPlus:
            return alpha_max_inverse_detonation_plus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseHybrid: 
            return alpha_max_inverse_hybrid(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDeflagration: 
            return alpha_max_inverse_deflagration(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::MixedHybridMinus:
            return alpha_max_mixed_hybrid_minus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::MixedHybridPlus:
            return alpha_max_mixed_hybrid_plus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
}

double alpha_min(
    HydroSolutionType sol_type, 
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    switch ( sol_type ) {
        case HydroSolutionType::Detonation:
            return alpha_min_detonation(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
       case HydroSolutionType::HybridPlus: 
            return alpha_min_hybrid_plus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
       case HydroSolutionType::Hybrid: 
            return alpha_min_hybrid(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::Deflagration: 
            return alpha_min_deflagration(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDetonation:
            return alpha_min_inverse_detonation(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDetonationPlus:
            return alpha_min_inverse_detonation_plus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseHybrid: 
            return alpha_min_inverse_hybrid(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDeflagration: 
            return alpha_min_inverse_deflagration(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::MixedHybridMinus:
            return alpha_min_mixed_hybrid_minus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::MixedHybridPlus:
            return alpha_min_mixed_hybrid_plus(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            );
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
}

double alpha_max(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( return_alpha_N ) {
        if ( xi_w <= cm ) {
            // direct deflagration
            return alpha_max_deflagration(
                xi_w, cp, cm, return_alpha_N, step_size, solve_mode
            ); 
        }
        // direct hybrid
        return alpha_max_hybrid(
            xi_w, cp, cm, return_alpha_N, step_size, solve_mode
        );
    }
    // if maximal alpha+ of inverse hybrid/deflagration (+ branch) below 1/3
    if ( compute_alpha_plus(cp, cm, cm) <= 1.0 / 3.0 ) { 
        // direct deflagration or direct hybrid ('-' branch)
        if ( xi_w < direct_jouguet_velocity(1.0 / 3.0, cm) ) { 
            return 1.0 / 3.0;
        }
        // direct detonation
        return compute_alpha_plus(xi_w, cm, cm);
    }
    // if maximal alpha+ of inverse hybrid/deflagration (+ branch) above 1/3
    // direct deflagration
    if ( xi_w <= get_v_minus_from_v_plus(cp, 1.0/3.0, RootBranch::Minus, cm) ) {
        return 1.0 / 3.0; // alpha(cp, xi_w) < 1/3
    }
    // inverse detonation on '+' branch or inverse hybrid
    if ( xi_w < cm ) {
        return compute_alpha_plus(cp, xi_w, cm);
    }
    // mixed hybrid
    if ( xi_w < cp ) {
        return compute_alpha_plus(cp, cm, cm);
    }
    // direct detonation
    return compute_alpha_plus(xi_w, cm, cm);
}

double alpha_min(
    double xi_w, 
    double cp,
    double cm,
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w < std::min(cm * cm, cp) ) {
        // inverse detonation
        return alpha_min_inverse_detonation(
            xi_w, cp, cm, return_alpha_N, step_size, solve_mode
        );
    } 
    if ( xi_w > cm * cm && xi_w < cp ) {
        // inverse hybrid
        return alpha_min_inverse_hybrid(
            xi_w, cp, cm, return_alpha_N, step_size, solve_mode
        );
    }
    if ( xi_w > std::max(cm * cm, cp) ) {
        // inverse deflagration
        return alpha_min_inverse_deflagration(
            xi_w, cp, cm, return_alpha_N, step_size, solve_mode
        );
    } 
    return alpha_min_deflagration(
        xi_w, cp, cm, return_alpha_N, step_size, solve_mode
    );
}

double inverse_detonation_plus_fold_critical_cm(double cp) {
    const double x = cp * cp;

    // Degree-20 fit to the critical line:
    // c_- = c_+ [1 - (1-c_+^2) P(c_+^2)].
    // RMSE in c_- = 2.720593257963e-04; max error = 5.383561285290e-04.
    static constexpr double coefficients[21] = {
        0.59554473655363,
        1.1148651488034,
        -81.2585134866,
        4126.3363162551,
        -109384.25963334,
        1755443.1006316,
        -18540650.356911,
        135837928.81916,
        -713851006.11834,
        2743787085.8081,
        -7766235707.9476,
        16044886879.26,
        -23297462068.555,
        20959991270.713,
        -4627060167.3229,
        -17097553634.584,
        28506898802.713,
        -23793821897.173,
        11881301806.434,
        -3386891784.8079,
        427063057.35081
    };

    double poly = coefficients[20];
    for (int k = 19; k >= 0; --k) {
        poly = poly * x + coefficients[k];
    }

    return cp * (1.0 - (1.0 - x) * poly);
}

bool check_noninjectivity(double cp, double cm) {
    return cp > cm && cm < inverse_detonation_plus_fold_critical_cm(cp);
}

HydroSolutionType classify_wave_fluid_frame(
    double xi_wall, 
    double alpha, 
    double cp,
    double cm,
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( ! check_velocity(xi_wall) ) { throw std::invalid_argument(
        "Invalid wall velocity: xi = " + std::to_string(xi_wall) + 
        ". Must be 0 < xi_wall < 1."
    ); }

    if (
        noninjectivity_check != CheckMode::SILENT &&
        check_noninjectivity(cp, cm) 
    ) {
        double cm_crit = inverse_detonation_plus_fold_critical_cm(cp);
        utils::check_failed(
            noninjectivity_check,
            "[CHECKS: non-injective hydrodynamic map] "
            "The chosen sound speeds c+ = " + std::to_string(cp)
            + ", c- = " + std::to_string(cm)
            + " (fitted fold-onset c- = " + std::to_string(cm_crit)
            + ") lie in the fitted region where the hydrodynamic map can "
            + "become non-injective for part of the parameter space. A fixed "
            + "(xi_w, alpha_N) may therefore admit multiple stationary "
            + "profiles. The automatic solver returns a single completion and "
            + "does not enumerate all branches or all roots of a folded "
            + "branch. Candidate branches can be constructed explicitly "
            + "through the branch-specific profile_calculator interface."
        );
    }

    // direct solutions 
    if ( alpha >= 0.0 ) {
        // direct deflagrations
        if ( xi_wall <= cm ) {
            if ( alpha < alpha_max_deflagration(
                xi_wall, cp, cm, true, step_size, solve_mode
            ) ) {
                return HydroSolutionType::Deflagration;
            }
        } else // xi_wall > cm
        // direct hybrids
        if ( xi_wall < direct_jouguet_velocity(alpha, cm) ) {   
            // '+' branch      
            if (
                cp > cm &&
                xi_wall < cp * cp / cm &&
                alpha < alpha_max_hybrid_plus(
                    xi_wall, cp, cm, true, step_size, solve_mode
                )
            ) {
                return HydroSolutionType::HybridPlus;
            }
            // '-' branch
            if ( alpha < alpha_max_hybrid(
                xi_wall, cp, cm, true, step_size, solve_mode
            ) ) {
                return HydroSolutionType::Hybrid;
            }
        } else // xi_wall >= xi_J > cm
        // direct detonations 
        if ( xi_wall >= cp ) {
            return HydroSolutionType::Detonation;
        } else // cm < xi_J <= xi_wall < cp => cp > cm
        // mixed hybrids
        if ( alpha >= alpha_min_mixed_hybrid_minus(
            xi_wall, cp, cm, true, step_size, solve_mode
        ) ) { 
            return HydroSolutionType::MixedHybridMinus; 
        } else
        // '+' branch
        if ( alpha >= alpha_min_mixed_hybrid_plus(
            xi_wall, cp, cm, true, step_size, solve_mode
        ) ) {
            return HydroSolutionType::MixedHybridPlus;
        } else
        // inverse hybrids reaching alpha_N >= 0 for c+ >> c-
        // (upper alpha_N bound of inverse hybrids is the same as the lower
        // bound of the mixed hybrid '+' branch)
        if (
            xi_wall > cm * cm &&
            alpha >= alpha_min_inverse_hybrid(
                xi_wall, cp, cm, true, step_size, solve_mode
            )
        ) {
            return HydroSolutionType::InverseHybrid;
        }
        // otherwise: no solution, fall through to the error at the end of the function
    // inverse solutions
    } else { // alpha < 0.0
        // inverse deflagrations
        if ( xi_wall >= std::max(cp, cm * cm) ) {
            if ( alpha > alpha_min_inverse_deflagration(
                xi_wall, cp, cm, true, step_size, solve_mode
            ) ) {
                return HydroSolutionType::InverseDeflagration;
            }
        } else // xi_wall < std::max(cp, cm*cm)
        if ( xi_wall <= std::min(cp, cm) ) {
            // inverse detonation, '-' branch
            if ( alpha > alpha_min_inverse_detonation(
                xi_wall, cp, cm, true, step_size, solve_mode
            ) ) {
                return HydroSolutionType::InverseDetonation;
            }
            if ( xi_wall >= cm * cm ) {
                // inverse detonation, '+' branch,
                if ( 
                    cp > cm &&
                    xi_wall >= xi_min(
                        HydroSolutionType::InverseDetonationPlus, cp, cm
                    ) &&
                    alpha > alpha_min_inverse_detonation_plus(
                        xi_wall, cp, cm, true, step_size, solve_mode
                    ) 
                ) {
                    return HydroSolutionType::InverseDetonationPlus;
                }
                // inverse hybrid
                if ( alpha > alpha_min_inverse_hybrid(
                    xi_wall, cp, cm, true, step_size, solve_mode
                ) ) {
                    return HydroSolutionType::InverseHybrid;
                }
            }
        } else // std::min(cp,cm) < xi_wall < std::max(cp, cm*cm) 
               //   => cm < xi_wall < cp || cm * cm > xi_wall > cp
        if (
            xi_wall < cp &&
            alpha >= alpha_min_inverse_hybrid(
                xi_wall, cp, cm, true, step_size, solve_mode
            )
        ) {
            // inverse hybrid
            if ( alpha <= alpha_min_mixed_hybrid_plus(
                xi_wall, cp, cm, true, step_size, solve_mode
            ) ) {
                return HydroSolutionType::InverseHybrid;
            }
            // mixed hybrid, v_+ = c_+
            if ( alpha <= alpha_min_mixed_hybrid_minus(
                xi_wall, cp, cm, true, step_size, solve_mode
            ) ) {
                return HydroSolutionType::MixedHybridPlus;
            }
            // mixed hybrid, v_- = c_-
            return HydroSolutionType::MixedHybridMinus;
        }
    }
    throw std::invalid_argument(
        "The parameters xi_wall = " + std::to_string(xi_wall) + 
        " and alpha_N = " + std::to_string(alpha) +
        " do not correspond to a valid hydrodynamic solution."
    );
}

HydroSolutionType classify_wave_wall_frame(
    double vp,
    double vm,
    double xi_w,
    double cp,
    double cm
) {
    if ( check_velocity(vp) && check_velocity(vm) && check_velocity(xi_w) ) {
        if (vp < vm) {
            if ( 
                vm < cm && 
                vm * vp < cp * cp &&
                are_equal(xi_w, vm)
            ) {
                return HydroSolutionType::Deflagration;
            }
            if ( 
                vp > cp && 
                vm * vp > cm * cm &&
                are_equal(xi_w, vp)
            ) {
                return HydroSolutionType::InverseDeflagration;
            }
            if (
                are_equal(vm, cm) &&
                xi_w > cm &&
                vp < xi_w &&
                xi_w * vp < cp * cp
            ) {
                return HydroSolutionType::Hybrid;
            }
        } else { // vp >= vm
            if ( 
                vm > cm && 
                vp > cp &&
                are_equal(xi_w, vp)
            ) {
                return HydroSolutionType::Detonation;
            }
            if ( 
                vm < cm && 
                vp < cp &&
                are_equal(xi_w, vm)
            ) {
                if (
                    cp > cm && vp > get_v_plus_extremum_from_v_minus(vm, cm)
                ) {
                    return HydroSolutionType::InverseDetonationPlus;
                }

                return HydroSolutionType::InverseDetonation;
            }
            if (
                are_equal(vm, cm) &&
                xi_w > cm &&
                xi_w > vp &&
                xi_w * vp < cp * cp
            ) {
                return HydroSolutionType::HybridPlus;
            }
            if (
                are_equal(vm, cm) &&
                xi_w > cm &&
                vp > xi_w &&
                vp < cp
            ) {
                return HydroSolutionType::MixedHybridMinus;
            }
            if (
                are_equal(vp, cp) &&
                xi_w < cp &&
                vm < xi_w &&
                vm > cm
            ) {
                return HydroSolutionType::MixedHybridPlus;
            }
        }
        if (
            are_equal(vp, cp) &&
            xi_w < cp &&
            vm > xi_w &&
            xi_w * vm > cm * cm
        ) {
            return HydroSolutionType::InverseHybrid;
        }
    }
    throw std::runtime_error(
        "Could not identify solution type. Invalid fluid velocities v+ = " 
        + std::to_string(vp) + " and v- = " + std::to_string(vm) + " with wall"
        " velocity xi = " + std::to_string(xi_w) + " and sound speeds c+ = "
        + std::to_string(cp) + " and c- = " + std::to_string(cm) + "."
    );
}

std::pair<double, double> compute_fluid_velocities(
    double xi_wall, 
    double alpha, 
    HydroSolutionType sol_type,
    double cp,
    double cm
) {
    double v_plus, v_minus;
    if ( debug ) { 
        std::cout << "[DEBUG: compute_fluid_velocities] (xi_w, alpha_+) = ("
                  << xi_wall << ", " << alpha << ") -> "; 
    }
    switch( sol_type ) {
        case HydroSolutionType::Deflagration: 
            v_minus = xi_wall;
            v_plus = get_v_plus_from_v_minus(
                v_minus, alpha, RootBranch::Minus, cm
            );
            break;
        case HydroSolutionType::Hybrid:
            v_minus = cm;
            v_plus = get_v_plus_from_v_minus(
                v_minus, alpha, RootBranch::Minus, cm
            );
            break;
        case HydroSolutionType::HybridPlus:
            v_minus = cm;
            v_plus = get_v_plus_from_v_minus(
                v_minus, alpha, RootBranch::Plus, cm
            );
            break;
        case HydroSolutionType::Detonation:
            v_plus = xi_wall; 
            v_minus = get_v_minus_from_v_plus(
                v_plus, alpha, RootBranch::Plus, cm
            );
            break;
        case HydroSolutionType::InverseDetonation:
            v_minus = xi_wall;
            v_plus = get_v_plus_from_v_minus(
                v_minus, alpha, RootBranch::Minus, cm
            );
            break;
        case HydroSolutionType::InverseDetonationPlus:
            v_minus = xi_wall;
            v_plus = get_v_plus_from_v_minus(
                v_minus, alpha, RootBranch::Plus, cm
            );
            break;
        case HydroSolutionType::InverseHybrid:
            v_plus = cp; 
            v_minus = get_v_minus_from_v_plus(
                v_plus, alpha, RootBranch::Plus, cm
            );
            break; 
        case HydroSolutionType::InverseDeflagration:
            v_plus = xi_wall; 
            v_minus = get_v_minus_from_v_plus(
                v_plus, alpha, RootBranch::Plus, cm
            );
            break; 
        case HydroSolutionType::MixedHybridMinus:
            v_minus = cm;
            v_plus = get_v_plus_from_v_minus(
                v_minus, alpha, RootBranch::Plus, cm
            );
            break;
        case HydroSolutionType::MixedHybridPlus:
            v_plus = cp;
            v_minus = get_v_minus_from_v_plus(
                v_plus, alpha, RootBranch::Plus, cm
            );
            break;
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
    if ( debug ) { 
        std::cout << sol_type << ": "
                  << "xi = " << xi_wall << ", alpha+ = " << alpha
                  << " c+ " << cp << ", c- = " << cm
                  << " -> v+ = " << v_plus << ", v- = " << v_minus 
                  << std::endl; 
    }
    // check if the velocities are valid 
    // (this should be fine for alpha_+ in the allowed range)
    if ( !(check_velocity(v_plus) && check_velocity(v_minus)) ) {
        throw std::runtime_error(
            "Invalid kinematic configuration: xi_w = " + std::to_string(xi_wall)
            + ", alpha_+ = " + std::to_string(alpha) + " (c+ = "
            + std::to_string(cp) + ", c- = " + std::to_string(cm) + "). Found "
            + "v_+ = " + std::to_string(v_plus) + " and v_- = "
            + std::to_string(v_minus) + "."
        );
    }
    // check if classification in terms of v_+ and v_- gives the correct result
    // note: for (inverse) hybrids, classification with xi_w instead of cs for 
    //       v_- (v_+) moves the solution into strong (inverse) deflagration 
    //       regime, but the code identifies these as hybrids
    if ( sol_type != classify_wave_wall_frame(
        v_plus, v_minus, xi_wall, cp, cm
    ) ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_fluid_velocities] "
            "Solution type changed after calculating the fluid velocities v+ = "
            + std::to_string(v_plus) + " and v- = " + std::to_string(v_minus)
            + " for xi_w = " + std::to_string(xi_wall) + " and alpha+ = " 
            + std::to_string(alpha)+ " (c+ = " + std::to_string(cp) + ", c- = "
            + std::to_string(cm) + "): " + to_string(sol_type) + " changed to " 
            + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, cp, cm
            )) + "."
        );
    }
    return {v_plus, v_minus};
}

double dxi_dv(double v, double xi, double cs) {
    double term1 = (xi - v);
    double term2 = (1.0 - xi * v);
    double term3 = (1.0 - v * v);
    double dxi_dt = xi * (term1 * term1 - cs * cs * term2 * term2);
    double dv_dt = 2.0 * v * cs * cs * term3 * term2;
    return dxi_dt / dv_dt;
}

double dv_dxi(double xi, double v, double cs) {
    double term1 = (xi - v);
    double term2 = (1.0 - xi * v);
    double term3 = (1.0 - v * v);
    double dxi_dt = xi * (term1 * term1 - cs * cs * term2 * term2);
    double dv_dt = 2.0 * v * cs * cs * term3 * term2;
    return  dv_dt / dxi_dt;
}

int integrate_xi_of_v(
    std::vector<double>& xi_vec, 
    std::vector<double>& v_vec, 
    double v_stop, 
    double dv, 
    bool stop_at_shock,
    double cs
) {
    if ( xi_vec.size() < 1 || xi_vec.size() != v_vec.size() ) {
        throw std::runtime_error(
            "The solution vectors must have equal size and contain the initial "
            "condition. Got sizes " + std::to_string(xi_vec.size()) + " and " 
            + std::to_string(v_vec.size()) + "."
        );
    }
    // use the last value stored in xi and v as initial conditions
    // (and append the results to the vectors)
    double xi = xi_vec.back(), v = v_vec.back();
    // decide whether to evolve in increasing or decreasing velocity direction
    if ( v_stop < v ) { dv = -std::abs(dv); } else { dv = std::abs(dv); }
    // determine when to stop (depending on the integration direction)
    std::function<bool(double)> advance_condition;
    if ( dv > 0.0 ) {
        advance_condition = [&](double v) { return v < v_stop; };
    } else {
        advance_condition = [&](double v) { return v > v_stop; };
    }
    if ( debug ) {
        std::cout << "[DEBUG: integrate_xi_of_v] Integrating "
                  << (dv > 0.0 ? "forward" : "backward") << " from (v, xi) = (" 
                  << v << ", " << xi << ") to v = " << v_stop << "." 
                  << std::endl;
    }
    double dxi, k1,k2,k3,k4;
    double shock1, shock2;
    // evolve differential equation with 4th-order Runge-Kutta
    while ( advance_condition(v + dv) ) {
        // determine step size in xi
        k1 = dv * dxi_dv(v, xi, cs);
        k2 = dv * dxi_dv(v + 0.5*dv, xi + 0.5*k1, cs);
        k3 = dv * dxi_dv(v + 0.5*dv, xi + 0.5*k2, cs);
        k4 = dv * dxi_dv(v + dv, xi + k3, cs);
        dxi = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        // stop if the self-similar coordinate becomes negative or super-luminal
        if ( !check_velocity(xi+dxi) ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_xi_of_v] Stopped at (v, xi) = ("
                          << v << ", " << xi << ") " << "as xi = " << xi+dxi 
                          << " becomes unphysical in next step." << std::endl;
            }
            return -1; 
        }
        // stop if the fluid velocity becomes super-luminal
        if ( std::abs(v+dv) >= 1.0 ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_xi_of_v] Stopped at (v, xi) = ("
                          << v << ", " << xi << ") " << "as v = " << v+dv 
                          << " becomes super-luminal in next step." 
                          << std::endl;
            }
            return -2; 
        }
        // stop if a shock front is crossed
        if ( stop_at_shock ) {
            shock1 = shock_condition(xi,v, cs);
            shock2 = shock_condition(xi+dxi, v+dv, cs);
            if ( shock1 * shock2 <= 0.0 ) {
                // crossed a shock: determine shock position by linear
                // interpolation and stop
                double t = shock1 / (shock1 - shock2);
                xi = xi + t * dxi;
                v = v + t * dv;
                xi_vec.emplace_back(xi);
                v_vec.emplace_back(v);
                if ( debug ) {
                    std::cout << "[DEBUG: integrate_xi_of_v] Encountered shock"
                              << " at (v, xi) = (" << v << ", " << xi << ")."
                              << std::endl;
                }
                return 1;
            }
        }
        // step forward
        xi += dxi; v += dv;
        xi_vec.emplace_back(xi);
        v_vec.emplace_back(v);
    }

    return 0;
}

int integrate_v_of_xi(
    std::vector<double>& xi_vec, 
    std::vector<double>& v_vec, 
    double v_stop, 
    double dxi, 
    bool stop_at_shock,
    double cs
) {
    if ( xi_vec.size() < 1 || xi_vec.size() != v_vec.size() ) {
        throw std::runtime_error(
            "The solution vectors must have equal size and contain the initial"
            " condition. Got sizes " + std::to_string(xi_vec.size()) + " and " 
            + std::to_string(v_vec.size()) + "."
        );
    }
    // use the last value stored in xi and v as initial conditions
    // (and append the results to the vectors)
    double xi = xi_vec.back(), v = v_vec.back();
    // decide whether to evolve in increasing or decreasing coordinate direction
    dxi = (v - v_stop) / dv_dxi(xi,v,cs) > 0.0 ? -std::abs(dxi) : std::abs(dxi);
    // determine when to stop (depending on the integration direction)
    std::function<bool(double)> advance_condition;
    if ( v < v_stop ) {
        advance_condition = [&](double v) { return v < v_stop; };
    } else {
        advance_condition = [&](double v) { return v > v_stop; };
    }
    if ( debug ) {
        std::cout << "[DEBUG: integrate_v_of_xi] Integrating " 
                  << (dxi > 0.0 ? "forward" : "backward") << " from (xi, v) = ("
                  << xi << ", " << v << ") to v = " << v_stop << "." 
                  << std::endl;
    }
    double dv, k1,k2,k3,k4;
    double shock1, shock2;
    double jouguet1 = jouguet_condition(xi, v, cs);
    double jouguet2;
    // evolve differential equation with 4th-order Runge-Kutta
    while ( advance_condition(v) ) {
        // determine step size in xi
        k1 = dxi * dv_dxi(xi, v, cs);
        k2 = dxi * dv_dxi(xi + 0.5*dxi, v + 0.5*k1, cs);
        k3 = dxi * dv_dxi(xi + 0.5*dxi, v + 0.5*k2, cs);
        k4 = dxi * dv_dxi(xi + dxi, v + k3, cs);
        dv = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        // stop if the self-similar coordinate becomes negative or super-luminal
        if ( !check_velocity(xi+dxi) ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_v_of_xi] Stopped at (xi, v) = ("
                          << xi << ", " << v << ") " << "as xi = " << xi+dxi 
                          << " becomes unphysical in next step." << std::endl;
            }
            return -1; 
        }
        // stop if the fluid velocity becomes super-luminal
        if ( std::abs(v+dv) >= 1.0 ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_v_of_xi] Stopped at (xi, v) = ("
                          << xi << ", " << v << ") " << "as v = " << v+dv 
                          << " becomes super-luminal in next step." 
                          << std::endl;
            }
            return -2; 
        }
        // stop if a shock front is crossed
        if ( stop_at_shock ) {
            shock1 = shock_condition(xi, v, cs);
            shock2 = shock_condition(xi+dxi, v+dv, cs);
            if ( shock1 * shock2 <= 0.0 ) {
                // crossed a shock: determine shock position by linear 
                // interpolation and stop, make sure that we don't overshoot
                // v_stop
                double t = (
                    advance_condition(v + dv) ?
                    shock1 / (shock1 - shock2) :
                    (v_stop - v)/dv
                );
                xi = xi + t * dxi;
                v = v + t * dv;
                xi_vec.emplace_back(xi);
                v_vec.emplace_back(v);
                if ( debug ) {
                    std::cout << "[DEBUG: integrate_v_of_xi] Encountered shock"
                              << " at (xi, vi) = (" << xi << ", " << v << ")." 
                              << std::endl;
                }
                return 1;
            }
        }
        // stop if the Jouguet front is crossed
        // this probably never occurs, the profile will probably evolve along 
        // the front, but anyway we are not going to solve the equations in this
        // direction.
        // jouguet1 = jouguet_condition(xi, v, cs);
        jouguet2 = jouguet_condition(xi+dxi, v+dv, cs);
        if ( jouguet1 * jouguet2 <= 0.0 ) {
            if ( debug ) {
                std::cout << "[DEBUG: integrate_v_of_xi] Stopped at (xi, v) = (" 
                          << xi << ", " << v << ")" << " as the Jouguet front"
                          << " is crossed in the next step." << std::endl;
            }
            return -3;
        }
        // step forward
        xi += dxi; v += dv;
        if ( !advance_condition(v) ) { break; }
        xi_vec.emplace_back(xi);
        v_vec.emplace_back(v);
    }
    return 0;
}

int integrate_velocity_profile(
    std::vector<double>& xi_vec,
    std::vector<double>& v_vec,
    double v_stop,
    double step_size,
    bool stop_at_shock,
    double cs,
    ProfileSolveMode solve_mode
) {
    switch ( solve_mode ) {
        case ProfileSolveMode::D_XI_D_V:
            return integrate_xi_of_v(
                xi_vec, v_vec, v_stop, step_size, stop_at_shock, cs
            );
        case ProfileSolveMode::D_V_D_XI:
            return integrate_v_of_xi(
                xi_vec, v_vec, v_stop, step_size, stop_at_shock, cs
            );
        case ProfileSolveMode::D_V_D_XI_FOR_SHOCKS:
            return (
                stop_at_shock ?
                integrate_v_of_xi(
                    xi_vec, v_vec, v_stop, step_size, stop_at_shock, cs
                ) :
                integrate_xi_of_v(
                    xi_vec, v_vec, v_stop, step_size, stop_at_shock, cs
                )
            );
    }
    throw std::runtime_error(
        "profile_calculator::integrate_velocity_profile got unknown solve mode."
    );
}

int compute_rarefaction_wave(
    std::vector<double>& xi,
    std::vector<double>& v,
    double c_s,
    double step_size,
    ProfileSolveMode solve_mode
) {
    double xi_wall = xi.back();
    double v_wall = v.back();
    double v_sonic = mu(xi_wall, c_s);
    if ( v_wall < std::min(0., v_sonic) || v_wall > std::max(0., v_sonic) ) {
        throw std::runtime_error(
            "Invalid initial conditions for rarefaction wave: (xi_w, v_w) = ("
            + std::to_string(xi_wall) + ", " + std::to_string(v_wall)
            + "), c_s = " + std::to_string(c_s) + " and mu(xi_w, c_s) = "
            + std::to_string(v_sonic) + "."
        );
    }
    double v_stop = xi_wall >= c_s ? v_zero : -v_zero;
    int status = integrate_velocity_profile(
        xi, v, v_stop, step_size, false, c_s, solve_mode
    );
    xi.push_back(c_s);
    v.push_back(0.);
    return status;
}

int compute_compression_wave(
    std::vector<double>& xi,
    std::vector<double>& v,
    double c_s,
    double step_size,
    ProfileSolveMode solve_mode
) {
    double xi_wall = xi.back();
    double v_wall = v.back();
    double v_sonic = mu(xi_wall, c_s);
    if ( v_wall > std::min(0., v_sonic) && v_wall < std::max(0., v_sonic) ) {
        throw std::runtime_error(
            "Invalid initial conditions for compression wave: (xi_w, v_w) = ("
            + std::to_string(xi_wall) + ", " + std::to_string(v_wall)
            + "), c_s = " + std::to_string(c_s) + " and mu(xi_w, c_s) = "
            + std::to_string(v_sonic) + "."
        );
    }
    double v_stop = v_wall > 0. ? v_zero : -v_zero;
    int status = integrate_velocity_profile(
        xi, v, v_stop, step_size, true, c_s, solve_mode
    );
    if (
        ( v_wall > 0.0 && xi.back() < c_s ) ||
        ( v_wall < 0.0 && xi.back() > c_s )
    ) {
        if ( debug ) {
            std::cout << "[DEBUG: compute_compression_wave] Integration ended"
                      << " (status " << status << ") at (xi, v) = ("
                      << xi.back() << ", " << v.back() << ") on the wrong side"
                      << " of the sound speed c_s = " << c_s << " for"
                      << " (xi_w, v_w) = (" << xi_wall << ", " << v_wall
                      << "). Appending (c_s, 0)." << std::endl;
        }
        xi.push_back(c_s); v.push_back(0.);
    }
    xi.push_back(xi.back());
    v.push_back(0.);
    return status;
}

size_t compute_velocity_profile_direct_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    double xi_wall = v_plus;
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus) 
        != HydroSolutionType::Detonation 
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and v- = "
            + std::to_string(v_minus) + " is not a detonation but of type " 
            + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0 = mu(xi_wall, v_minus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Direct Detonation:"
                  << " xi_w = " << xi_wall << ", v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to sound speed (behind the wall)
    status = compute_rarefaction_wave(xi, v, c_minus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_direct_detonation] "
            "Expected integration return status 0 (reached v_stop), but got "
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add first element
    xi.push_back(0.);
    v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add last two elements
    xi.push_back(xi.back()); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_direct_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    double xi_wall = v_minus;
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus) 
        != HydroSolutionType::Deflagration 
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and v- = " 
            + std::to_string(v_minus) + " is not a deflagration but of type " 
            + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0 = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Direct Deflagration:"
                  << " xi_w = " << xi_wall << ", v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add first two elements
    xi.push_back(0.); xi.push_back(xi_wall);
    v.push_back(0.); v.push_back(0.);
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to shock (in front of the wall)
    status = compute_compression_wave(xi, v, c_plus, step_size, solve_mode);
    if ( status != 1 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_direct_deflagration] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last element
    xi.push_back(1.);
    v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_direct_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double xi_wall, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    double v_minus = c_minus;
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus)
        != HydroSolutionType::Hybrid 
        &&
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus)
        != HydroSolutionType::HybridPlus
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and xi_w = " 
            + std::to_string(xi_wall) + " is not a hybrid but of type " 
            + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }

    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0m = mu(xi_wall, v_minus);
    double v0p = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Direct Hybrid: xi_w = " 
                  << xi_wall << ", v(xi_w^-) = " << v0m << ", v(xi_w^+) = " 
                  << v0p << "." << std::endl;
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0m);
    // integrate from wall to sound speed (behind the wall)
    status = compute_rarefaction_wave(xi, v, c_minus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_direct_hybrid] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" + 
            (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", xi_w = " 
            + std::to_string(xi_wall) + "."
        );
    }
    // add first element
    xi.push_back(0.);
    v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0p);
    // integrate from wall to shock (in front of the wall)
    status = compute_compression_wave(xi, v, c_plus, step_size, solve_mode);
    if ( status != 1 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_direct_hybrid] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", xi_w = " 
            + std::to_string(xi_wall) + "."
        );
    }
    // add last element
    xi.push_back(1.);
    v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_inverse_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus,
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    double xi_wall = v_minus;
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus) 
        != HydroSolutionType::InverseDetonation 
        &&
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus) 
        != HydroSolutionType::InverseDetonationPlus
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v- = "
            + std::to_string(v_minus) + " is not an inverse detonation but of "
            "type " + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0 = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Inverse Detonation:"
                  << " xi_w = " << xi_wall << " and v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add first two elements
    xi.push_back(0.); xi.push_back(xi_wall);
    v.push_back(0.); v.push_back(0.);
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to sound speed (in front of the wall)
    status = compute_rarefaction_wave(xi, v, c_plus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_inverse_detonation] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last element
    xi.push_back(1.);
    v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_inverse_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    double xi_wall = v_plus;
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus) != 
        HydroSolutionType::InverseDeflagration 
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and v- = " 
            + std::to_string(v_minus) + " is not an inverse deflagration but of"
            " type " + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }

    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0 = mu(xi_wall, v_minus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Inverse Deflagration:"
                  << " xi_w = " << xi_wall << " and v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to shock (behind the wall)
    status = compute_compression_wave(xi, v, c_minus, step_size, solve_mode);
    if ( status != 1 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_inverse_deflagration] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add first element
    xi.push_back(0.);
    v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add last two elements
    xi.push_back(xi.back()); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_inverse_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double v_minus, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    double v_plus = c_plus;
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus) 
        != HydroSolutionType::InverseHybrid 
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and v- = "
            + std::to_string(v_minus) + " is not an inverse hybrid but of type "
            + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0m = mu(xi_wall, v_minus);
    double v0p = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Inverse Hybrid: xi_w = "
                  << xi_wall << ", v(xi_w^-) = " << v0m << ", v(xi_w^+) = " 
                  << v0p << "." << std::endl;
    }

    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0m);
    // integrate from wall to shock (behind the wall)
    status = compute_compression_wave(xi, v, c_minus, step_size, solve_mode);
    if ( status != 1 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_inverse_hybrid] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for xi_w = " + std::to_string(xi_wall) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add first element
    xi.push_back(0.);
    v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0p);
    // integrate from wall to sound speed (in front of the wall)
    status = compute_rarefaction_wave(xi, v, c_plus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_inverse_hybrid] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for xi_w = " + std::to_string(xi_wall) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last element
    xi.push_back(1.);
    v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_mixed_hybrid_minus(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double xi_wall, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    double v_minus = c_minus;
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus)
        != HydroSolutionType::MixedHybridMinus 
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and xi_w = " 
            + std::to_string(xi_wall) + " is not a mixed hybrid with v- = c-"
            " but of type " + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }

    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0m = mu(xi_wall, v_minus);
    double v0p = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Mixed Hybrid (v- = c-):"
                  << " xi_w = " << xi_wall << ", v(xi_w^-) = " << v0m
                  << ", v(xi_w^+) = " << v0p << "." << std::endl;
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0m);
    // integrate from wall to sound speed (behind the wall)
    status = compute_rarefaction_wave(xi, v, c_minus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_mixed_hybrid_minus] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" + 
            (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", xi_w = " 
            + std::to_string(xi_wall) + "."
        );
    }
    // add first element
    xi.push_back(0.);
    v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0p);
    // integrate from wall to speed of sound (in front of the wall)
    status = compute_rarefaction_wave(xi, v, c_plus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_mixed_hybrid_minus] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for xi_w = " + std::to_string(xi_wall) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last element
    xi.push_back(1.);
    v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_mixed_hybrid_plus(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double v_minus, 
    double c_plus,
    double c_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    double v_plus = c_plus;
    if ( 
        classify_wave_wall_frame(v_plus, v_minus, xi_wall, c_plus, c_minus)
        != HydroSolutionType::MixedHybridPlus 
    ) {
        throw std::runtime_error(
            "The solution with v+ = " + std::to_string(v_plus) + " and xi_w = " 
            + std::to_string(xi_wall) + " is not a mixed hybrid with v+ = c+"
            " but of type " + to_string(classify_wave_wall_frame(
                v_plus, v_minus, xi_wall, c_plus, c_minus
            )) + " for c+ = " + std::to_string(c_plus) + " and c- = " 
            + std::to_string(c_minus) + "."
        );
    }

    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0m = mu(xi_wall, v_minus);
    double v0p = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Mixed Hybrid (v+ = c+):"
                  << " xi_w = " << xi_wall << ", v(xi_w^-) = " << v0m 
                  << ", v(xi_w^+) = " << v0p << "." << std::endl;
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0m);
    // integrate from wall to sound speed (behind the wall)
    status = compute_rarefaction_wave(xi, v, c_minus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_mixed_hybrid_plus] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" + 
            (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v+ = " + std::to_string(v_plus) + ", xi_w = " 
            + std::to_string(xi_wall) + "."
        );
    }
    // add first element
    xi.push_back(0.);
    v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0p);
    // integrate from wall to speed of sound (in front of the wall)
    status = compute_rarefaction_wave(xi, v, c_plus, step_size, solve_mode);
    if ( status != 0 ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile_mixed_hybrid_plus] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for xi_w = " + std::to_string(xi_wall) + ", v- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last element
    xi.push_back(1.);
    v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t get_shock_index(HydroSolutionType sol_type, size_t size) {
    switch (sol_type) {
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
        case HydroSolutionType::HybridPlus:
            return size - 2;
        case HydroSolutionType::InverseDeflagration:
        case HydroSolutionType::InverseHybrid:
            return 2;
        default:
            break;
    }
    throw std::runtime_error(
        "Solution type " + to_string(sol_type) + " does not have a shock front."
    );
}

std::pair<double, double> get_shock_front(
    HydroSolutionType sol_type, 
    const std::vector<double>& xi, 
    const std::vector<double>& v
) {
    check_profiles(xi, v);
    size_t idx = get_shock_index(sol_type, xi.size());
    if ( 
        sol_type == HydroSolutionType::Deflagration 
        || sol_type == HydroSolutionType::Hybrid 
        || sol_type == HydroSolutionType::HybridPlus
    ) { --idx; } // for direct transition, return v_-, else v_+
    return {xi[idx], v[idx]};
}

size_t compute_velocity_profile(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_N,
    double cp,
    double cm,
    double step_size,
    ProfileSolveMode solve_mode
) {
    HydroSolutionType sol_type = classify_wave_fluid_frame(
        xi_wall, alpha_N, cp, cm, step_size, solve_mode
    );
    double alpha_plus = find_alpha_plus(
        sol_type, alpha_N, xi_wall, cp, cm, step_size, solve_mode
    );
    return compute_velocity_profile(
        sol_type, xi, v, xi_wall, alpha_plus, cp, cm, step_size, solve_mode
    );
}

size_t compute_velocity_profile(
    HydroSolutionType sol_type, 
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_p, 
    double cp,
    double cm,
    double step_size,
    ProfileSolveMode solve_mode
) {
    
    // solve matching conditions at bubble wall
    auto [v_plus, v_minus] = compute_fluid_velocities(
        xi_wall, alpha_p, sol_type, cp, cm
    );

    size_t idx_wall; // index of first element outside the bubble
    switch (sol_type) {
        case HydroSolutionType::Deflagration:
            idx_wall = compute_velocity_profile_direct_deflagration(
                xi, v, v_plus, v_minus, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::Detonation:
            idx_wall = compute_velocity_profile_direct_detonation(
                xi, v, v_plus, v_minus, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::Hybrid:
        case HydroSolutionType::HybridPlus:
            idx_wall = compute_velocity_profile_direct_hybrid(
                xi, v, v_plus, xi_wall, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::InverseDeflagration:
            idx_wall = compute_velocity_profile_inverse_deflagration(
                xi, v, v_plus, v_minus, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::InverseDetonation:
        case HydroSolutionType::InverseDetonationPlus:
            idx_wall = compute_velocity_profile_inverse_detonation(
                xi, v, v_plus, v_minus, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::InverseHybrid:
            idx_wall = compute_velocity_profile_inverse_hybrid(
                xi, v, xi_wall, v_minus, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::MixedHybridMinus:
            idx_wall = compute_velocity_profile_mixed_hybrid_minus(
                xi, v, v_plus, xi_wall, cp, cm, step_size, solve_mode
            );
            break;
        case HydroSolutionType::MixedHybridPlus:
            idx_wall = compute_velocity_profile_mixed_hybrid_plus(
                xi, v, xi_wall, v_minus, cp, cm, step_size, solve_mode
            );
            break;
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) + "."
            );
    }

    if ( !are_equal(xi_wall, xi[idx_wall]) ) {
        utils::check_failed(
            consistency_check,
            "[CHECKS: compute_velocity_profile] "
            "The wall velocity of the calculated velocity spectrum does not"
            " match the input value: " + std::to_string(xi[idx_wall]) + " vs. " 
            + std::to_string(xi_wall) + "."
        );
    }

    // check shock front
    double v_shock = v[2];      // for inverse deflagration/hybrid
    double dv = v_shock - v[3]; // for inverse deflagration/hybrid
    switch (sol_type) {
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
        case HydroSolutionType::HybridPlus:
            v_shock = v[v.size()-3];
            dv = v[v.size()-4] - v_shock;
        case HydroSolutionType::InverseDeflagration:
        case HydroSolutionType::InverseHybrid:
            if ( std::abs(v_shock) < std::abs(dv) ) {
                utils::check_failed(
                    consistency_check,
                    "[CHECKS: compute_velocity_profile] v_shock " 
                    + std::to_string(v_shock) 
                    + " is too small or zero. Increase precision!"
                );
            }
            break;
        case HydroSolutionType::Detonation:
        case HydroSolutionType::InverseDetonation:
        case HydroSolutionType::InverseDetonationPlus:
        case HydroSolutionType::MixedHybridMinus:
        case HydroSolutionType::MixedHybridPlus:
        default:
            break;
    }

    return idx_wall;
}

// TODO: can be split into constant and linear part. Constant part needs to be 
//       evaluated at start/end and jumps only.
double integrate_enthalpy_segment(
    double xi1,
    double xi2,
    double v1,
    double v2,
    double cs
) {
    if ( xi1 == xi2 ) { return 0.0; } // jump position
    if ( v1 == v2 ) { return 0.0; } // integral vanishes if v is constant
    // approximate xi(v) as linear function: xi(v) = xi0 + dxi * v
    double xi0 = (xi1 * v2 - xi2 * v1) / (v2 - v1);
    double dxi = (xi2 - xi1) / (v2 - v1);
    // partial fraction decomposition:
    // gamma^2 * mu = (xi0-(1-dxi)*v) / (1-v^2) / (1-xi0*v-dxi*v^2) 
    //              = - v/(1-v^2) - (dxi*v + xi0)/(dxi*v^2+xi0*v-1))
    // first integral:
    double result = 0.5 * (std::log(1.0 - v2*v2) - std::log(1.0 - v1*v1));
    // second integral:
    result += 0.5 * (std::log(1.0 - xi1*v1) - std::log(1.0 - xi2*v2));
    double delta = 4.0 * dxi + xi0*xi0;
    if ( delta > 0 ) { 
        result += 0.5 * xi0 / std::sqrt(delta) * std::log(
            (2.0 * dxi * v1 + xi0 - std::sqrt(delta)) 
            / (2.0 * dxi * v1 + xi0 + std::sqrt(delta)) 
            * (2.0 * dxi * v2 + xi0 + std::sqrt(delta))
            / (2.0 * dxi * v2 + xi0 - std::sqrt(delta))
        );
    } else {
        result += xi0 / std::sqrt(-delta) * (
            std::atan((2.0 * dxi * v1 + xi0) / std::sqrt(-delta))
            - std::atan((2.0 * dxi * v2 + xi0) / std::sqrt(-delta))
        );
    }
    return (1.0 + 1.0 / cs / cs) * result;
}

std::vector<double> compute_enthalpy_density(
    const std::vector<double>& xi, 
    const std::vector<double>& v, 
    size_t idx_wall, 
    size_t idx_shock,
    double cp,
    double cm
) {
    check_profiles(xi, v);
    const size_t N = xi.size();
    std::vector<double> w(N);
    w[N-1] = 1.0;
    double cs = cp; // start in front of the wall
    for ( size_t i=xi.size()-1; i>0; --i ) {
        // note: if idx_shock=-1, its actual value is PTRDIFF_MAX-1
        //       (since size_t is unsigned) but this is fine here
        if ( i == idx_wall || i == idx_shock ) { 
            if ( xi[i] != xi[i-1] ) { 
                throw std::runtime_error(
                    std::to_string(i) + " is not a jump position: No jump"
                    " between (xi,v) = (" + std::to_string(xi[i-1]) + ", " 
                    + std::to_string(v[i-1]) + ") and (xi, v) = ()"
                    + std::to_string(xi[i]) + ", " + std::to_string(v[i]) + "."
                );
            }
            // jump in enthalpy density: go to rest frame of the front and match
            double vp = mu(xi[i], v[i]);
            double vm = mu(xi[i], v[i-1]);
            if ( !check_velocity(vp) || !check_velocity(vm) ) {
                throw std::runtime_error(
                    "Encountered invalid velocities after boosting back to wall"
                    " frame: v_+ = " + std::to_string(vp) + " and v_- = " 
                    + std::to_string(vm) + "."
                );
            }
            if ( debug ) { 
                std::cout << "[DEBUG: compute_enthalpy_density] reconstructed"
                          << " (v_+, v-)  = (" << vp << ", " << vm << ") at "
                          << (i==idx_wall ? "bubble" : "shock") << " front." 
                          << std::endl;
            }
            // apply enthalpy jump across then front
            w[i-1] = w[i] * vp / vm * (1.0 - vm * vm) / (1.0 - vp * vp);
            // change cs if we crossed the bubble wall
            if ( i == idx_wall ) { cs = cm; }
        } else {
            w[i-1] = w[i] * std::exp(integrate_enthalpy_segment(
                xi[i], xi[i-1], v[i], v[i-1], cs
            ));
        }
    }
    return w;
}

std::vector<double> compute_energy_density_profile(
    const std::vector<double>& w,
    size_t idx_wall,
    double epsilon,
    double cp,
    double cm
) {
    std::vector<double> e(w.size());
    double r = 1.0 / (1.0 + cm * cm );
    double epsp = 0.0, epsm = 0.0;
    if ( epsilon > 0.0 ) { epsp = epsilon; } else { epsm = -epsilon; }
    for ( size_t i=0; i<idx_wall; ++i )  { 
        e[i] = r * w[i] + epsm; 
    }
    r = 1.0 / (1.0 + cp * cp );
    for ( size_t i=idx_wall; i<w.size(); ++i )  { 
        e[i] = r * w[i] + epsp; 
    }
    return e;
}

std::vector<double> compute_pressure_profile(
    const std::vector<double>& w,
    size_t idx_wall,
    double epsilon,
    double cp,
    double cm
) {
    std::vector<double> p(w.size());
    double r = cm * cm / (1.0 + cm * cm );
    double epsp = 0.0, epsm = 0.0;
    if ( epsilon > 0.0 ) { epsp = epsilon; } else { epsm = -epsilon; }
    for ( size_t i=0; i<idx_wall; ++i )  { 
        p[i] = r * w[i] - epsm; 
    }
    r = cp * cp / (1.0 + cp * cp );
    for ( size_t i=idx_wall; i<w.size(); ++i )  { 
        p[i] = r * w[i] - epsp; 
    }
    return p;
}

std::vector<double> compute_temperature_profile(
    const std::vector<double>& w,
    size_t idx_wall,
    double T_N,
    double Psi_N,
    double cp,
    double cm
) {
    std::vector<double> T(w.size());
    double nu_p = get_nu_from_cs(cp);
    double nu_m = get_nu_from_cs(cm); 
    double w_N = w[w.size()-1]; 
    for ( size_t i=0; i<idx_wall; ++i )  { 
        T[i] = T_N * std::pow(w[i] / w_N / Psi_N, 1.0 / nu_m); 
    }
    for ( size_t i=idx_wall; i<w.size(); ++i )  { 
        T[i] = T_N * std::pow(w[i] / w_N, 1.0 / nu_p); 
    }
    return T;
}

std::vector<double> compute_entropy_density_profile(
    const std::vector<double>& w, 
    const std::vector<double>& T
) {
    check_profiles(T, w);
    std::vector<double> s(w.size());
    for ( size_t i=0; i<w.size(); ++i )  { 
        s[i] = w[i] / T[i];
    }
    return s;
}

size_t check_entropy_flux(
    const std::vector<double>& xi,
    const std::vector<double>& v,
    const std::vector<double>& s
) {
    check_profiles(xi, v);
    check_profiles(xi, s);
    double v_minus, v_plus, j_minus, j_plus;
    for ( size_t i = 0; i < xi.size()-1; ++i ) {
        if ( std::abs(xi[i+1] - xi[i]) > jump_tolerance ) {
            continue;
        }
        v_minus = mu(xi[i], v[i]);
        v_plus = mu(xi[i], v[i+1]);
        j_minus = s[i] * v_minus * gamma(v_minus);
        j_plus = s[i+1] * v_plus * gamma(v_plus);
        if ( j_minus < j_plus && !are_equal(j_minus, j_plus) ) {
            if ( debug ) {
                std::cout << "[DEBUG: check_entropy_flux] Profile violates "
                          << "entropy flux condition for discontinuity at xi = "
                          << xi[i+1] << " (index " << i+1 << "): " << j_minus
                          << " < " << j_plus << "." << std::endl;
            }
            return i+1;
        }
    }
    return NO_INDEX;
}

double compute_alpha_bag_from_pseudotrace(double alpha, double cp, double cm) {
    double cp2 = cp * cp;
    double cm2 = cm * cm;
    double alpha0 = compute_alpha_plus(cp2, cm2, cm);
    return (alpha - alpha0) * 3.0 * cm2 * (1.0 + cp2) / (1.0 + cm2);
}

double compute_alpha_pseudotrace_from_bag(double alpha, double cp, double cm) {
    double cp2 = cp * cp;
    double cm2 = cm * cm;
    double alpha0 = compute_alpha_plus(cp2, cm2, cm);
    return alpha0 + alpha * (1.0 + cm2) / 3.0 / cm2 / (1.0 + cp2);
}

double rescale_alpha_by_enthalpy(
    double alpha,
    double r_w,
    double cp,
    double cm
) {
    double alpha0 = compute_alpha_plus(cp * cp, cm * cm, cm);
    return alpha0 * (1.0 - r_w) + r_w * alpha;
}

double compute_epsilon(double w, double alpha, double cp, double cm) {
    double alpha_bag = compute_alpha_bag_from_pseudotrace(alpha, cp, cm);
    return alpha_bag * w / (1.0 + cp * cp);
}

double compute_alpha_N(
    HydroSolutionType sol_type, 
    double alpha_plus, 
    double xi_wall, 
    double cp,
    double cm,
    double step_size,
    ProfileSolveMode solve_mode
) {
    auto [v_plus, v_minus] = compute_fluid_velocities(
        xi_wall, alpha_plus, sol_type, cp, cm
    ); 
    std::vector<double> xi{xi_wall}, v{mu(xi_wall, v_plus)};
    double xish, mush;
    double matching_factor = 1.0;
    switch (sol_type) {
        // for detonation and inverse deflagration:
        //     alpha_N and alpha_+ are identical
        case HydroSolutionType::Detonation:
        case HydroSolutionType::InverseDeflagration:
            return alpha_plus;
        // for deflagration and hybrids: 
        //     evolve w from bubble wall to shock front
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
        case HydroSolutionType::HybridPlus:
            integrate_velocity_profile(
                xi, v, v_zero, step_size, true, cp, solve_mode
            );
            xish = xi.back(); // vsh+
            mush = mu(xish, v.back()); // vsh-
            if ( !check_velocity(mush) ) {
                throw std::runtime_error(
                    "Encountered invalid velocity after boosting back to shock"
                    " frame: mu(xi_sh, v_sh) = " + std::to_string(mush) + "."
                );
            }
            matching_factor =  (
                xish / mush * (1.0 - mush * mush) / (1.0 - xish*xish)
            ); // matching at shock front
            break;
        // for inverse detonation and inverse hybrids: evolve w from bubble wall
        // to speed-of-sound
        case HydroSolutionType::InverseHybrid:
        case HydroSolutionType::InverseDetonation:
        case HydroSolutionType::InverseDetonationPlus:
        case HydroSolutionType::MixedHybridPlus:
        case HydroSolutionType::MixedHybridMinus:
            integrate_velocity_profile(
                xi, v, -v_zero, step_size, false, cp, solve_mode
            );
            break;
    }
    double log_w_N_over_w_plus = 0.0;
    for ( size_t i=0; i<xi.size()-1; ++i ) { 
        log_w_N_over_w_plus += integrate_enthalpy_segment(
            xi[i+1], xi[i], v[i+1], v[i], cp
        ); 
    }
    return rescale_alpha_by_enthalpy(
        alpha_plus, 
        matching_factor * std::exp(log_w_N_over_w_plus), 
        cp, 
        cm
    );
}

double find_alpha_plus(
    double alpha_N, 
    double xi_wall, 
    double cp,
    double cm,
    double step_size, 
    ProfileSolveMode solve_mode
) {
    HydroSolutionType sol_type = classify_wave_fluid_frame(
        xi_wall, alpha_N, cp, cm, step_size, solve_mode
    ); 
    return find_alpha_plus(
        sol_type, alpha_N, xi_wall, cp, cm, step_size, solve_mode
    );
}

double find_alpha_plus(
    HydroSolutionType sol_type,
    double alpha_N, 
    double xi_wall, 
    double cp,
    double cm,
    double step_size, 
    ProfileSolveMode solve_mode
) {    
    // check if xi_wall and alphaN are valid
    double alpha_N_min = alpha_min(
        sol_type, xi_wall, cp, cm, true, step_size, solve_mode
    );
    double alpha_N_max = alpha_max(
        sol_type, xi_wall, cp, cm, true, step_size, solve_mode
    );
    if ( alpha_N < alpha_N_min || alpha_N > alpha_N_max ) { 
        throw std::runtime_error(
            "Invalid combination of wall velocity and transition strength:"
            " xi_w = " + std::to_string(xi_wall) + ", alpha_N = " 
            + std::to_string(alpha_N) + " for " + to_string(sol_type)
            + " solutions (c+ = " + std::to_string(cp) + ", c- = " 
            + std::to_string(cm) + ")."
        ); 
    }
    // bracket alpha+ based on solution type
    double alpha_p_min = alpha_min(sol_type, xi_wall, cp, cm, false);
    double alpha_p_max = alpha_max(sol_type, xi_wall, cp, cm, false);
    alpha_p_min += alpha_abs_tol;
    alpha_p_max -= alpha_abs_tol;
    double alpha_fix = compute_alpha_plus(cp * cp, cm * cm , cm);
    switch ( sol_type ) {
        // note: for direct detonation and inverse deflagration, alphaN and
        //       alpha+ are the same
        case HydroSolutionType::Detonation:
        case HydroSolutionType::InverseDeflagration:
            return alpha_N; 
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
        case HydroSolutionType::HybridPlus:
            if ( alpha_N < alpha_fix ) { // alpha+ > alphaN
                alpha_p_min = std::max(alpha_p_min, alpha_N);
            } else { // alpha+ < alphaN 
                alpha_p_max = std::min(alpha_p_max, alpha_N);
            }
            break;
        case HydroSolutionType::InverseDetonation:
        case HydroSolutionType::InverseDetonationPlus:
        case HydroSolutionType::InverseHybrid:
        case HydroSolutionType::MixedHybridMinus:
        case HydroSolutionType::MixedHybridPlus:
            if ( alpha_N > alpha_fix ) { // alpha+ > alphaN
                alpha_p_min = std::max(alpha_p_min, alpha_N);
            } else { // alpha+ < alphaN 
                alpha_p_max = std::min(alpha_p_max, alpha_N);
            }
            break;
    }
    
    if ( debug ) {
        std::cout << "[DEBUG: find_alpha_plus] " << sol_type << ": bracketed"
                  << " alpha in [" << alpha_p_min << ", " << alpha_p_max 
                  << "] for" << " (xi_w, alpha_N) = " << xi_wall << ", " 
                  << alpha_N << ")." << std::endl;
    }

    // find root of the difference (absolute, so that alpha_N = 0 works)
    auto alpha_abs_difference = [
        sol_type, alpha_N, xi_wall, cp, cm, step_size, solve_mode
    ](double alpha_plus) {
        return compute_alpha_N(
            sol_type, alpha_plus, xi_wall, cp, cm, step_size, solve_mode
        ) - alpha_N;
    };
    // use combination of absolute and relative tolerance
    auto tolerance = [
        rel_tol = alpha_rel_tol, abs_tol = alpha_abs_tol
    ](double a, double b) {
        return std::abs(b - a) <= std::max(
            abs_tol, 
            rel_tol*std::max(std::abs(a), std::abs(b))
        );
    };

    // check that alpha_N really changes sign in the bracketing interval
    if ( 
        alpha_abs_difference(alpha_p_max) * 
        alpha_abs_difference(alpha_p_min) > 0. 
    ) {
        throw std::runtime_error(
            "Cannot find alpha_+: No sign change in initial bracketing range ("
            + std::to_string(alpha_p_min) + ", " + std::to_string(alpha_p_max) 
            + "): alpha_N = (" + std::to_string(compute_alpha_N(
                sol_type, alpha_p_min, xi_wall, cp, cm, step_size, solve_mode
            )) + ", " + std::to_string(compute_alpha_N(
                sol_type, alpha_p_max, xi_wall, cp, cm, step_size, solve_mode
            )) + ")."
        );
    }

    // find root
    std::uintmax_t N_iter = 1000;
    auto result = boost::math::tools::toms748_solve(
        alpha_abs_difference, alpha_p_min, alpha_p_max, tolerance, N_iter
    );
    double alpha_plus = .5*(result.first+result.second);

    if ( debug ) {
        double delta_alpha_plus = alpha_abs_difference(alpha_plus);
        std::cout << "[DEBUG: find_alpha_plus] Found alpha_+ = " << alpha_plus 
                  << ". " << "Absolute difference in alpha_N: " 
                  << delta_alpha_plus << "." << std::endl;
    }

    return alpha_plus;
}

double compute_efficiency(
    const std::vector<double>& xi,
    const std::vector<double>& v,
    const std::vector<double>& w,
    size_t idx_wall,
    double alpha_N,
    double cp,
    double cm, 
    EfficiencyNormalization normalization
) {
    size_t size = xi.size();
    check_profiles(xi, v); 
    check_profiles(xi, w);
    std::vector<double> integrand(size);
    double xi_wall = xi[idx_wall];
    double xi_wall3 = xi_wall * xi_wall * xi_wall;
    double w_N = w.back();
    double epsilon = compute_epsilon(w_N, alpha_N, cp, cm);
    // find size of the perturbed region
    size_t i_max = size - 1;
    while ( i_max > 0 && v[i_max] == 0.0 ) { --i_max; }
    double xi_max3 = xi[i_max] * xi[i_max] * xi[i_max];
    // set prefactor
    double prefactor = 0.0;
    switch (normalization) {
        case EfficiencyNormalization::BubbleVolume:
            prefactor = 3.0 / xi_wall3;
            break;
        case EfficiencyNormalization::TotalEnergyDensity:
            prefactor = 3.0 / ( 
                xi_max3 * w_N/(1.0 + cp*cp) 
                + xi_wall3 * std::max(0.0, epsilon)
            );
            break;
        case EfficiencyNormalization::VacuumEnergy:
            prefactor = (
                epsilon > 0.0 ?
                3.0 / xi_wall3 / epsilon :
                3.0 * (1.0 + cp*cp) / w_N / xi_max3
            );
            break;
        case EfficiencyNormalization::AlphaPseudoTrace:
            prefactor = (
                alpha_N > 0.0 ?
                4.0 / (alpha_N * w_N * xi_wall3 ):
                3.0 * (1.0 + cp*cp) / w_N / xi_max3
            );
            break; 
    }
    if ( prefactor == 0.0 ) { 
        throw std::runtime_error(
            "Bad efficiency normalization: vanishing prefactor."
        );
    }
    // integrate kinetic energy
    #pragma omp simd
    for ( size_t i = 0; i < size; ++i ) {
        double temp = v[i] * gamma(v[i]) * xi[i];
        integrand[i] = w[i] * temp * temp;
    }
    return prefactor * utils::integrate_trapez(xi, integrand);
}

double compute_enthalpy_ratio_LTE(
    double xi_wall, 
    double alpha_N, 
    double cp,
    double cm,
    double step_size, 
    ProfileSolveMode solve_mode
) {
    HydroSolutionType sol_type = classify_wave_fluid_frame(
        xi_wall, alpha_N, cp, cm, step_size, solve_mode
    );
    double alpha_plus = find_alpha_plus(
        sol_type, alpha_N, xi_wall, cp, cm, step_size, solve_mode
    );
    auto [v_plus, v_minus] = compute_fluid_velocities(
        xi_wall, alpha_plus, sol_type, cp, cm
    );
    double nu = get_nu_from_cs(cm);
    return v_plus / v_minus * std::pow(
        (1.0 - v_plus * v_plus) / (1.0 - v_minus * v_minus), 0.5 * nu - 1.0
    );
}

// TODO: CHECK CAREFULLY!
double find_LTE_wall_velocity(
    double alpha_N, 
    double Psi, 
    double xi_min, 
    double xi_max, 
    double cp,
    double cm,
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_min >= xi_max ) { throw std::runtime_error(
        "Cannot find LTE wall velocity: xi_min and xi_max are not properly"
        " ordered: xi_min = " + std::to_string(xi_min) + " > xi_max = " 
        + std::to_string(xi_max) + "."
    ); }
    if ( 
        alpha_N > 0 
        && xi_min < direct_jouguet_velocity(alpha_N, cm) 
        && xi_max > direct_jouguet_velocity(alpha_N, cm) 
    ) {
        throw std::runtime_error(
            "Cannot find LTE wall velocity: The bracketing interval (" 
            + std::to_string(xi_min) + ", " + std::to_string(xi_max) 
            + ") crosses the Jouguet velocity " 
            + std::to_string(direct_jouguet_velocity(alpha_N, cm)) + "."
        );
    } else if ( 
        alpha_N < 0 
        && xi_min < std::max(cp, cm * cm) 
        && xi_max > std::max(cp, cm * cm) 
    ) {
        throw std::runtime_error(
            "Cannot find LTE wall velocity: The bracketing interval (" 
            + std::to_string(xi_min) + ", " + std::to_string(xi_max) 
            + ") crosses boundary of inverse deflagration solutions."
        );
    }
    if ( debug ) {
        std::cout << "[DEBUG: find_LTE_wall_velocity] Bracketing LTE wall"
                  << " velocity in [" << xi_min << ", " << xi_max << "] for"
                  << " (alpha_N, Psi) = (" << alpha_N << ", " << Psi << ")."
                  << std::endl;
    }
    // find root of the relative difference
    auto psi_rel_difference = [alpha_N, Psi, cp, cm, step_size, solve_mode](
        double xi_wall
    ) { 
        return compute_enthalpy_ratio_LTE(
            xi_wall, alpha_N, cp, cm, step_size, solve_mode
        ) / Psi - 1.0; 
    };
    // use combination of absolute and relative tolerance
    auto tolerance = [rel_tol = xi_wall_rel_tol, abs_tol = xi_wall_abs_tol](
        double a, double b
    ) {
        return std::abs(b - a) <= std::max(
            abs_tol, 
            rel_tol*std::max(std::abs(a), std::abs(b))
        );
    };
    // check that Psi really changes sign in the bracketing interval
    if ( psi_rel_difference(xi_min) * psi_rel_difference(xi_max) > 0. ) {
        throw std::runtime_error(
            "Cannot find LTE wall velocity: No sign change in initial"
            " bracketing range (" + std::to_string(xi_min) + ", " 
            + std::to_string(xi_max) + "): Delta_Psi/Psi = ("
            + std::to_string(psi_rel_difference(xi_min)) + ", " 
            + std::to_string(psi_rel_difference(xi_max)) + ")."
        );
    }

    // find root
    std::uintmax_t N_iter = 1000;
    auto result = boost::math::tools::toms748_solve(
        psi_rel_difference, xi_min, xi_max, tolerance, N_iter
    );
    double xi_wall = .5*(result.first+result.second);

    if ( debug ) {
        double delta_Psi = psi_rel_difference(xi_wall);
        std::cout << "[DEBUG: find_LTE_wall_velocity] Found xi_wall = " 
                  << xi_wall << ". " << "Relative difference in Psi: " 
                  << delta_Psi << "." << std::endl;
    }

    return xi_wall;
}

bool check_profiles(
    const std::vector<double>& xi, 
    const std::vector<double>& profile, 
    bool throw_error
) {
    size_t size = xi.size();
    if ( size != profile.size() || size < 6) {
        std::string message = 
            "The self-similar coordinate and the profile vectors must have the"
            " same size and at least size 6. Got sizes " + std::to_string(size) 
            + " and " + std::to_string(profile.size()) + ".";

        if ( throw_error ) { throw std::runtime_error(message); }
        if ( debug ) { 
            std::cerr << "[DEBUG: check_profiles] " << message << std::endl; 
        }
        return false;
    }
    return true;
}

} // namespace profile_calculator
