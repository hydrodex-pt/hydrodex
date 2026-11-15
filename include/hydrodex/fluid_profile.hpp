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

#include "hydrodex/constants.hpp"
#include "hydrodex/profile_calculator.hpp"

#include <vector>

namespace hydrodex {

/**
 * @brief Class representing the fluid profile of a bubble in a cosmological
 *     phase transition.
 */
class FluidProfile {
public:

    /// @name Constructors
    ///@{

    /**
     * @brief Constructor for the \c FluidProfile from precomputed profiles.
     * @param xi_ Self-similar coordinate \f$\xi = r/t\f$
     * @param v_ Fluid velocity profile \f$v(\xi)\f$
     * @param w_ Enthalpy density profile \f$w(\xi)\f$
     * @param e_ Energy density profile \f$e(\xi)\f$
     * @param p_ Pressure profile \f$p(\xi)\f$
     * @param idx_w Index of the first element outside the bubble
     * @param idx_sh Index of the first element outside the shock front
     *     (use \c constants::NO_INDEX if there is no shock front, default:
     *     \c constants::NO_INDEX)
     * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
     *     \c constants::cs_bag)
     * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
     *     \c constants::cs_bag)
     */
    FluidProfile(
        const std::vector<double>& xi_, 
        const std::vector<double>& v_, 
        const std::vector<double>& w_, 
        const std::vector<double>& e_, 
        const std::vector<double>& p_, 
        size_t idx_w, 
        size_t idx_sh = constants::NO_INDEX,
        double cp = constants::cs_bag,
        double cm = constants::cs_bag
    );

    /** 
     * @brief Constructor for the \c FluidProfile from the wall velocity
     *     \f$\xi_w\f$ and the pseudotrace phase transition strength parameter
     *     \f$\alpha_N\f$.
     * @param xi_w Wall velocity \f$\xi_w\f$
     * @param alpha_N Pseudotrace strength parameter \f$\alpha_N\f$
     *     (far in front of the wall)
     * @param cp Speed of sound \f$c_+\f$ ahead of the wall (default:
     *     \c constants::cs_bag)
     * @param cm Speed of sound \f$c_-\f$ behind the wall (default:
     *     \c constants::cs_bag)
     * @param step_size Step size for the fluid velocity (or coordinate) when
     *     solving the differential equation (default:
     *     \c profile_calculator::profile_step_size)
     * @param solve_mode \c ProfileSolveMode switch determining which
     *     differential equation for the velocity profile is solved (default:
     *     \c profile_calculator::profile_solve_mode)
     */
    FluidProfile(
        double xi_w, 
        double alpha_N, 
        double cp = constants::cs_bag,
        double cm = constants::cs_bag,
        double step_size=profile_calculator::profile_step_size,
        profile_calculator::ProfileSolveMode solve_mode 
            = profile_calculator::profile_solve_mode
    );

    ///@}

    /** 
     * @brief Get the self-similar coordinate \f$\xi=r/t\f$ of the profile.
     * @return Reference to the vector containing the coordinate
     */
    const std::vector<double>& get_xi_coordinate() const;

    /** 
     * @brief Get the fluid velocity profile \f$v(\xi)\f$.
     * @return Reference to the vector containing the velocity
     */
    const std::vector<double>& get_fluid_velocity() const;

    /** 
     * @brief Get the enthalpy density profile \f$w(\xi)\f$.
     * @return Reference to the vector containing the enthalpy
     */
    const std::vector<double>& get_enthalpy_density() const;

    /** 
     * @brief Get the energy density profile \f$e(\xi)\f$.
     * @return Reference to the vector containing the energy
     */
    const std::vector<double>& get_energy_density() const;

    /** 
     * @brief Get the pressure profile \f$p(\xi)\f$.
     * @return Reference to the vector containing the pressure
     */
    const std::vector<double>& get_pressure() const;


    /** 
     * @brief Get the temperature profile \f$T(\xi)\f$.
     * @param Psi_N Enthalpy ratio \f$\Psi_N\f$ across the wall at nucleation
     * @param T_N Nucleation temperature \f$T_N\f$ (default: 1.0)
     * @return Vector containing the temperature profile
     */
    std::vector<double> get_temperature(
        double Psi_N,
        double T_N = 1.0
    ) const;
    
    /** 
     * @brief Get the entropy density profile \f$s(\xi)\f$.
     * @param Psi_N Enthalpy ratio \f$\Psi_N\f$ across the wall at nucleation
     * @param s_N Entropy density \f$s_N\f$ at nucleation (default: 1.0)
     * @return Vector containing the entropy density profile
     */
    std::vector<double> get_entropy_density(
        double Psi_N,
        double s_N = 1.0
    ) const;

    /**
     * @brief Get the index of the element in \f$\xi\f$ corresponding to the
     *     bubble wall position.
     * @return Index of the bubble wall position (first element outside the
     *     wall)
     */
    size_t get_wall_index() const;

    /**
     * @brief Get the index of the element in \f$\xi\f$ corresponding to the
     *     shock front position.
     * @return Index of the shock front position (first element outside the
     *     shock), or \c constants::NO_INDEX if there is no shock front
     */
    size_t get_shock_index() const;

    /**
     * @brief Get the wall velocity \f$\xi_w\f$.
     * @return Value of the wall velocity \f$\xi_w\f$
     */
    double get_wall_velocity() const;

    /**
     * @brief Get the pseudotrace phase transition strength \f$\alpha_N\f$.
     * @return Value of the transition strength \f$\alpha_N\f$
     */
    double get_transition_strength() const;

    /**
     * @brief Get the speed of sound \f$c_+\f$ in the false vacuum phase.
     * @return Speed of sound \f$c_+\f$ in the false vacuum phase
     */
    double get_sound_speed_plus() const;

    /**
     * @brief Get the speed of sound \f$c_-\f$ in the true vacuum phase.
     * @return Speed of sound \f$c_-\f$ in the true vacuum phase
     */
    double get_sound_speed_minus() const;

    /**
     * @brief Get the wall-frame fluid velocity \f$v_+\f$ ahead of the wall.
     * @return Wall-frame fluid velocity \f$v_+\f$ ahead of the wall
     */
    double get_v_plus() const;

    /**
     * @brief Get the wall-frame fluid velocity \f$v_-\f$ behind the wall.
     * @return Wall-frame fluid velocity \f$v_-\f$ behind the wall
     */
    double get_v_minus() const;

    /**
     * @brief Get the pseudotrace phase transition strength \f$\alpha_+\f$ ahead
     *     of the wall.
     * @return Pseudotrace phase transition strength \f$\alpha_+\f$ ahead of the
     *     bubble wall
     */
    double get_wall_transition_strength() const;

    /**
     * @brief Get the self-similar coordinate \f$\xi_\text{sh}\f$ of the shock
     *     front.
     * @return Position \f$\xi_\text{sh}\f$ of the shock front, or `NaN` if
     *     there is no shock front
     */
    double get_shock_position() const;

    /**
     * @brief Get the type of the hydrodynamic solution.
     * @return Hydrodynamic solution type
     */
    profile_calculator::HydroSolutionType get_solution_type() const;

    /**
     * @brief Check whether the profile has a shock front. 
     * @return \c true if the profile has a shock front, \c false otherwise
     */
    bool has_shock() const;

    /**
     * @brief Get the mean enthalpy density \f$\bar{w}\f$ far in front of
     *     the wall.
     * @return Mean enthalpy density \f$\bar{w}\f$
     */
    double get_mean_enthalpy_density() const;

    /** 
     * @brief Get the efficiency \f$K\f$ (or \f$\kappa_\mathrm{sw}\f$).
     * @param normalization \c EfficiencyNormalization switch to determine how
     *     to normalize the efficiency (default:
     *     \c EfficiencyNormalization::TotalEnergyDensity)
     * @return Efficiency \f$K\f$ or \f$\kappa_\mathrm{sw}\f$
     */
    double get_efficiency(
        profile_calculator::EfficiencyNormalization normalization
            = profile_calculator::EfficiencyNormalization::TotalEnergyDensity
    ) const;

    /**
     * @brief Set verbosity.
     * @param b Boolean value to which the verbosity flag is set
     */
    static void set_verbose(bool b);

    /**
     * @brief Get verbosity.
     * @return Boolean value to which the verbosity flag is set
     */
    static bool is_verbose();

private:
    inline static bool verbose = false; ///< verbosity flag
    std::vector<double> xi; ///< self-similar coordinate \f$\xi=r/t\f$
    std::vector<double> v;  ///< fluid velocity profile \f$v(\xi)\f$
    std::vector<double> w;  ///< enthalpy density profile \f$w(\xi)\f$
    std::vector<double> e;  ///< energy density profile \f$e(\xi)\f$
    std::vector<double> p;  ///< pressure profile \f$p(\xi)\f$
    size_t idx_wall;        ///< index of the first element outside the wall
    size_t idx_shock;       ///< index of the first element outside the shock
    double xi_wall;         ///< wall velocity \f$\xi_w\f$
    double alpha;           ///< pseudotrace transition strength \f$\alpha_N\f$
    double c_plus;          ///< speed of sound \f$c_+\f$ outside the wall
    double c_minus;         ///< speed of sound \f$c_-\f$ inside the wall
    /// type of the hydrodynamic solution
    profile_calculator::HydroSolutionType sol_type; 

    /**
     * @brief Check consistency and validity of the fluid profile data.
     * @throws std::runtime_error if the check fails
     */
    void check_profiles() const;

}; // class FluidProfile

} // namespace hydrodex