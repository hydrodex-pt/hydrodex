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

#include "hydrodex/settings.hpp"

#include <gsl/gsl_errno.h>

namespace hydrodex::settings {
    double v_zero = 0.0;
    double jump_tolerance = 1e-9;
    bool debug = false;
    CheckMode consistency_check = CheckMode::WARNING;
    CheckMode noninjectivity_check = CheckMode::WARNING;

    void set_gsl_abort_on_error(bool abort) {
        if ( abort ) { gsl_set_error_handler(nullptr); } // GSL default
        else { gsl_set_error_handler_off(); }
    }
} // namespace settings