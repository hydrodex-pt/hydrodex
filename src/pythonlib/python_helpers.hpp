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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>
#include <stdexcept>

#include "hydrodex/constants.hpp"

using namespace hydrodex;

/** 
 * @brief Convert numpy array to std::vector<double>.
 * @param arr 1D numpy array to convert
 * @returns Converted array as std::vector<double>
 */
inline std::vector<double> np2vec(const pybind11::array_t<double>& arr) {
    if (arr.ndim() != 1) {
        throw std::runtime_error("Input array must be 1D.");
    }
    // unchecked access respects the strides, so non-contiguous arrays (e.g.
    // slices like x[::2]) are copied correctly as well
    auto view = arr.unchecked<1>();
    std::vector<double> vec(view.shape(0));
    for (pybind11::ssize_t i = 0; i < view.shape(0); ++i) { vec[i] = view(i); }
    return vec;
}

/**
 * @brief Convert std::vector<double> to numpy array.
 * @param vec std::vector<double> to convert
 * @returns Converted vector as 1D numpy array
 */
inline pybind11::array_t<double> vec2np(const std::vector<double>& vec) {
    return pybind11::array_t<double>(vec.size(), vec.data());
}

/**
 * @brief Recover the dynamic default value for an argument.
 * @details Converts the input \p obj to the same data type as the default value
 *     \p default_value or returns \p default_value if \p obj is \c None.
 * @param obj Python object holding the argument value (or \c None)
 * @param default_value Default value used if \p obj is \c None
 * @returns \p obj (as the correct data type) or \p default_value
 */
template <typename T>
inline T recover_dynamic_default(
    const pybind11::object& obj,
    const T& default_value
) {
    return obj.is_none() ? default_value : obj.cast<T>();
}

/**
 * @brief Convert Python object to an index (\c size_t).
 * @details Converts the index (integer) to \c size_t and \c None to
 *     \c constants::NO_INDEX.
 * @param obj Python object (integer) holding the index
 * @returns Index as \c size_t, mapping \c None to \c constants::NO_INDEX
 */
inline size_t obj2idx(pybind11::object obj) {
    return ( 
        obj.is_none() ? constants::NO_INDEX : obj.cast<size_t>() 
    );
}

/**
 * @brief Convert \c size_t to integer index in Python.
 * @details Converts the \c size_t index to an integer index in Python and
 *     \c constants::NO_INDEX to \c None.
 * @param idx \c size_t index in C++
 * @returns Integer index in Python (or \c None if the @p idx is 
 *     \c constants::NO_INDEX)
 */
inline pybind11::object idx2obj(size_t idx) {
    if ( idx == constants::NO_INDEX ) { return pybind11::none(); }
    return pybind11::int_(idx);
}

