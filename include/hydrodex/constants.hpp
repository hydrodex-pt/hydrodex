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

#include <cstddef>

/**
 * @namespace hydrodex::constants
 * @brief Numerical constants used throughout HydróDex.
 */
namespace hydrodex::constants {
    /// \f$\pi\f$
    inline constexpr double Pi = 3.141592653589793;
    /// speed of sound in the MIT bag model, \f$c_s = \frac{1}{\sqrt{3}}\f$
    inline constexpr double cs_bag = 0.57735026918962576451;       
    /// index value used when not applicable
    /// (e.g. when no shock front is present)
    inline constexpr size_t NO_INDEX = static_cast<size_t>(-1); 
} // namespace constants