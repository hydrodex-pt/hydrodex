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

#include "hydrodex/fluid_profile.hpp"

#include "hydrodex/constants.hpp"
#include "hydrodex/profile_calculator.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace hydrodex {

FluidProfile::FluidProfile(
    const std::vector<double>& xi_, 
    const std::vector<double>& v_, 
    const std::vector<double>& w_, 
    const std::vector<double>& e_, 
    const std::vector<double>& p_,
    size_t idx_w, 
    size_t idx_sh,
    double cp,
    double cm
) : 
    xi(xi_), 
    v(v_), 
    w(w_), 
    e(e_), 
    p(p_), 
    idx_wall(idx_w), 
    idx_shock(idx_sh), 
    xi_wall(xi_[idx_w]),
    c_plus(cp),
    c_minus(cm)
{ 
    if ( verbose ) { 
        std::cout << "Initializing fluid profile from data. " << std::endl; 
    }
    // check profiles
    check_profiles(); 
    if ( verbose ) { std::cout << "Profiles are okay." << std::endl; }
    // classify solution
    double v_plus = get_v_plus();
    double v_minus = get_v_minus();
    sol_type = profile_calculator::classify_wave_wall_frame(
        v_plus, v_minus, xi_wall, c_plus, c_minus
    );
    double alpha_plus = get_wall_transition_strength();
    double alpha_N = profile_calculator::rescale_alpha_by_enthalpy(
        alpha_plus, w[idx_wall] / w.back(), c_plus, c_minus
    );
    alpha = alpha_N;
    if ( verbose ) {
        std::cout << "Classified profile as " << sol_type << " solution with"
                  << " v_+ = " << v_plus << " and v_- = " << v_minus << ","
                  << " alpha_+ = " << alpha_plus << ".\nThis corresponds to "
                  << "xi_w = " << xi_wall << " and alpha_N = " << alpha_N << "."
                  << std::endl; 
    }
    if ( ( idx_shock == constants::NO_INDEX ) != ( 
        sol_type == profile_calculator::HydroSolutionType::Detonation || 
        sol_type == profile_calculator::HydroSolutionType::MixedHybridMinus ||
        sol_type == profile_calculator::HydroSolutionType::MixedHybridPlus ||
        sol_type == profile_calculator::HydroSolutionType::InverseDetonation ||
        sol_type == profile_calculator::HydroSolutionType::InverseDetonationPlus
    ) ) { throw std::runtime_error(
        "The solution is of type " + to_string(sol_type) + " but " + (
            idx_shock == constants::NO_INDEX
            ? "no shock index" : "shock index " + std::to_string(idx_shock)
        ) + " was provided."
    ); } 
}

FluidProfile::FluidProfile(
    double xi_w, 
    double alpha_N, 
    double cp,
    double cm,
    double step_size, 
    profile_calculator::ProfileSolveMode solve_mode
) : xi_wall(xi_w), alpha(alpha_N), c_plus(cp), c_minus(cm) {
    if ( verbose ) {
        std::cout << "Calculating fluid profile for xi_w = " << xi_w 
                  << " and alpha_N = " << alpha_N << "." << std::endl;
    }

    sol_type = profile_calculator::classify_wave_fluid_frame(
        xi_wall, alpha_N, c_plus, c_minus, step_size, solve_mode
    );
    double alpha_plus = profile_calculator::find_alpha_plus(
        sol_type, alpha_N, xi_wall, c_plus, c_minus, step_size, solve_mode
    );
    auto [vplus, vminus] = profile_calculator::compute_fluid_velocities(
        xi_wall, alpha_plus, sol_type, c_plus, c_minus
    );

    if ( verbose ) {
        std::cout << "Found " << sol_type << " solution with alpha_+ = " 
                  << alpha_plus << ".\n" << "Fluid velocities: v_+ = " << vplus
                  << ", v_- = " << vminus << ".\nCalculating velocity profile."
                  << std::endl;
    }

    idx_wall = profile_calculator::compute_velocity_profile(
        sol_type, 
        xi, 
        v, 
        xi_wall, 
        alpha_plus, 
        c_plus, 
        c_minus, 
        step_size, 
        solve_mode
    );
    idx_shock = constants::NO_INDEX;
    if ( 
        sol_type != profile_calculator::HydroSolutionType::Detonation &&
        sol_type != profile_calculator::HydroSolutionType::MixedHybridMinus &&
        sol_type != profile_calculator::HydroSolutionType::MixedHybridPlus &&
        sol_type != profile_calculator::HydroSolutionType::InverseDetonation &&
        sol_type != profile_calculator::HydroSolutionType::InverseDetonationPlus
    ) {
        idx_shock = profile_calculator::get_shock_index(sol_type, xi.size());
    }
    if ( verbose ) {
        std::cout << "Velocity profile calculated with " << xi.size() 
                  << " points.\n  Velocity jumps from " << v[idx_wall-1] 
                  << " to " << v[idx_wall] << " at the bubble wall, xi = " 
                  << xi[idx_wall] << " (index: " << idx_wall << ")." ;
        if ( idx_shock != constants::NO_INDEX ) { 
            std::cout << "\n  Velocity jumps from " << v[idx_shock-1] << " to " 
                      << v[idx_shock] << " at the shock position, xi = " 
                      << xi[idx_shock] << " (index: " << idx_shock << ").";
        }
        std::cout << std::endl;
    }
    w = profile_calculator::compute_enthalpy_density(
        xi, v, idx_wall, idx_shock, c_plus, c_minus
    );
    double wN = w.back();
    if (verbose) { 
        std::cout << "Enthalpy density calculated. wN: " << wN << std::endl; 
    }

    double epsilon = profile_calculator::compute_epsilon(
        wN, alpha_N, c_plus, c_minus
    );
    if ( verbose ) { 
        std::cout << "Vacuum energy: epsilon = " << epsilon << std::endl; 
    }
    e = profile_calculator::compute_energy_density_profile(
        w, idx_wall, epsilon, c_plus, c_minus
    );
    p = profile_calculator::compute_pressure_profile(
        w, idx_wall, epsilon, c_plus, c_minus
    );

    if ( verbose ) { std::cout << "All profiles calculated." << std::endl; }
}

const std::vector<double>& FluidProfile::get_xi_coordinate() const 
    { return xi; }
const std::vector<double>& FluidProfile::get_fluid_velocity() const 
    { return v; }
const std::vector<double>& FluidProfile::get_enthalpy_density() const 
    { return w; }
const std::vector<double>& FluidProfile::get_energy_density() const 
    { return e; }
const std::vector<double>& FluidProfile::get_pressure() const 
    { return p; }

std::vector<double> FluidProfile::get_temperature(
    double Psi_N, double T_N
) const {
    return profile_calculator::compute_temperature_profile(
        w, idx_wall, T_N, Psi_N, c_plus, c_minus
    );
}
std::vector<double> FluidProfile::get_entropy_density(
    double Psi_N, double s_N
) const {
    double T_N = w[w.size()-1] / s_N;
    std::vector<double> T = profile_calculator::compute_temperature_profile(
        w, idx_wall, T_N, Psi_N, c_plus, c_minus
    );
    return profile_calculator::compute_entropy_density_profile(w, T);
}


size_t FluidProfile::get_wall_index() const { return idx_wall; }
size_t FluidProfile::get_shock_index() const { return idx_shock; }
double FluidProfile::get_wall_velocity() const { return xi_wall; }
double FluidProfile::get_transition_strength() const { return alpha; }
double FluidProfile::get_sound_speed_plus() const { return c_plus; }
double FluidProfile::get_sound_speed_minus() const { return c_minus; }
double FluidProfile::get_v_plus() const {
    return profile_calculator::mu(xi_wall, v[idx_wall]);
}
double FluidProfile::get_v_minus() const {
    return profile_calculator::mu(xi_wall, v[idx_wall-1]);
}
double FluidProfile::get_wall_transition_strength() const {
    double v_plus = get_v_plus();
    double v_minus = get_v_minus();
    return profile_calculator::compute_alpha_plus(v_plus, v_minus, c_minus);
}
double FluidProfile::get_shock_position() const {
    if ( idx_shock == constants::NO_INDEX ) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return xi[idx_shock];
}
profile_calculator::HydroSolutionType FluidProfile::get_solution_type() const
    { return sol_type; }
bool FluidProfile::has_shock() const 
    { return idx_shock != constants::NO_INDEX; }
double FluidProfile::get_mean_enthalpy_density() const { return w.back(); }

double FluidProfile::get_efficiency(
    profile_calculator::EfficiencyNormalization normalization
) const {
    return profile_calculator::compute_efficiency(
        xi, v, w, idx_wall, alpha, c_plus, c_minus, normalization
    );
}

void FluidProfile::set_verbose(bool b) { verbose = b; };
bool FluidProfile::is_verbose() { return verbose; };

void FluidProfile::check_profiles() const {
    profile_calculator::check_profiles(xi, v);
    profile_calculator::check_profiles(xi, w);
    profile_calculator::check_profiles(xi, e);
    profile_calculator::check_profiles(xi, p);
    if ( idx_wall < 2 || idx_wall > xi.size() - 2 ) { throw std::runtime_error(
        "Invalid value for the wall velocity index: " + std::to_string(idx_wall)
        + ". The index must be between 2 and " + std::to_string(xi.size()-2)
        + " (for size " + std::to_string(xi.size()) + ")." 
    );}
    if ( ! profile_calculator::are_equal(xi[idx_wall-1], xi[idx_wall]) ) { 
        throw std::runtime_error(
            "Wall index " + std::to_string(idx_wall) + " does not correspond to"
            " the bubble wall. " + "The self-similar coordinate changes from " 
            + std::to_string(xi[idx_wall-1]) + " to " 
            + std::to_string(xi[idx_wall]) + "."
        );
    }
    if ( idx_shock != constants::NO_INDEX ) {
        if ( idx_shock < 2 || idx_shock > xi.size() - 2 ) { 
            throw std::runtime_error(
                "Invalid value for the shock front index: " 
                + std::to_string(idx_shock) + ". The index must be between 2"
                " and " + std::to_string(xi.size()-2) + " (for size " 
                + std::to_string(xi.size()) + ")." 
            );
        }
        if ( ! profile_calculator::are_equal(xi[idx_shock-1], xi[idx_shock]) ) {
            throw std::runtime_error(
                "Shock index " + std::to_string(idx_shock) + " does not"
                " correspond to the shock front. The self-similar coordinate"
                " changes from " + std::to_string(xi[idx_shock-1]) + " to " + 
                std::to_string(xi[idx_shock]) + "."
            );
        }  
    }
}

} // namespace hydrodex