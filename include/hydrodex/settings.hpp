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

/**
 * @namespace hydrodex::settings
 * @brief Global settings controlling numerical tolerances, debug output and
 *     the handling of failed checks in HydróDex.
 */
namespace hydrodex::settings {

    /// report mode options for non-critical consistency checks
    enum class CheckMode {
        ERROR,   ///< throw exception if a check fails
        WARNING, ///< show a warning if a check fails
        SILENT   ///< ignore failed checks
    };

    /// small number corresponding to zero velocity when calculating profiles
    extern double v_zero;
    /// tolerance for considering two values of the coordinate as equal
    extern double jump_tolerance;
    /// debug output flag
    extern bool debug;
    /// report mode for non-critical consistency checks
    extern CheckMode consistency_check;
    /// report mode for non-injectivity checks
    extern CheckMode noninjectivity_check;

    /**
     * @brief Switch GSL's default error handler (abort on error) on or off.
     *
     * With @p abort = false, failed GSL integrations are reported as
     * exceptions instead of terminating the program. 
     *
     * @param abort Whether GSL should abort the program on errors
     * 
     * @note This is a global GSL setting. It also affects other code using the
     *     same GSL library.
     */
    void set_gsl_abort_on_error(bool abort);

} // namespace settings