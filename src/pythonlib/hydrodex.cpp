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

#include <gsl/gsl_errno.h>

#include "hydrodex/version.hpp"
#include "hydrodex_bindings.hpp"

namespace py = pybind11;

// =============================================
// HydróDex module
// =============================================
PYBIND11_MODULE(hydrodex, m) {
    m.doc() = "HydróDex\n\n"
        "A C++ code for the calculation of the fluid profiles for direct and"
        " inverse cosmological phase transitions and the corresponding"
        " gravitational wave signal within the sound shell model.\n\n\n"
        "References\n"
        "----------\n"
        "  [Barni+ (2027)]\n"
        "    G. Barni and E. Madge,\n"
        "    \"HydroDex: gotta match 'em all!\"\n"
        "    [arXiv:2611.17411 [physics.flu-dyn]].\n"
        "  [Barni+ (2026)]\n"
        "    G. Barni, S. Blasi, E. Madge and M. Vanvlasselaer,\n"
        "    \"Gravitational waves from the sound shell model: direct and"
        " inverse phase transitions in the early Universe,\"\n"
        "    JCAP 03 (2026), 066 [arXiv:2510.21439 [hep-ph]].\n"
        "  [Barni+ (2025)]\n"
        "    G. Barni, S. Blasi and M. Vanvlasselaer,\n"
        "    \"Inverse bubbles from broken symmetries,\"\n"
        "    Phys. Rev. D 112 (2025) 9, 095006 [arXiv:2503.01951 [hep-ph]].\n"
        "  [Giombi+ (2025)]\n"
        "    L. Giombi, J. Dahl and M. Hindmarsh,\n"
        "    \"Signatures of the speed of sound on the gravitational wave power"
        " spectrum from sound waves,\"\n"
        "    JCAP 01 (2025), 100 [arXiv:2409.01426 [gr-qc]].\n"
        "  [Barni+ (2024)]\n"
        "    G. Barni, S. Blasi and M. Vanvlasselaer,\n"
        "    \"The hydrodynamics of inverse phase transitions,\"\n"
        "    JCAP 10 (2024), 042 [arXiv:2406.01596 [hep-ph]].\n"
        "  [Roper+ (2024)]\n"
        "    A. Roper Pol, S. Procacci and C. Caprini,\n"
        "    \"Characterization of the gravitational wave spectrum from sound"
        " waves within the sound shell model,\"\n"
        "    Phys. Rev. D 109 (2024) 6, 063531 [arXiv:2308.12943 [gr-qc]].\n"
        "  [Wang+ (2022)]\n"
        "    X. Wang, F. P. Huang and Y. Li,\n"
        "    \"Sound velocity effects on the phase transition gravitational"
        " wave spectrum in the sound shell model,\"\n"
        "    Phys. Rev. D 105 (2022) 10, 103513 [arXiv:2112.14650"
        " [astro-ph.CO]].\n"
        "  [Guo+ (2021)]\n"
        "    H.-K. Guo, K. Sinha, D. Vagie and G. White,\n"
        "    \"Phase Transitions in an Expanding Universe: Stochastic"
        " Gravitational Waves in Standard and Non-Standard Histories,\"\n"
        "    JCAP 01 (2021), 001 [arXiv:2007.08537 [hep-ph]].\n"
        "  [Giese+ (2020)]\n"
        "    F. Giese, T. Konstandin and J. van de Vis,\n"
        "    \"Model-independent energy budget of cosmological first-order"
        " phase transitions - A sound argument to go beyond the bag model,\"\n"
        "    JCAP 07 (2020), 057 [arXiv:2004.06995 [astro-ph.CO]].\n"
        "  [Hindmarsh+ (2019)]\n"
        "    M. Hindmarsh and M. Hijazi,\n"
        "    \"Gravitational waves from first order cosmological phase"
        " transitions in the Sound Shell Model,\"\n"
        "    JCAP 12 (2019), 062 [arXiv:1909.10040 [astro-ph.CO]].\n"
        "  [Hindmarsh (2018)]\n" 
        "    M. Hindmarsh,\n"
        "    \"Sound shell model for acoustic gravitational wave production at"
        " a first-order phase transition in the early Universe,\"\n"
        "    Phys. Rev. Lett. 120 (2018) no. 7, 071301 [arXiv:1608.04735"
        " [astro-ph.CO]].\n"
        "  [Leitao+ (2015)]\n"
        "    L. Leitao and A. Megevand,\n"
        "    \"Hydrodynamics of phase transition fronts and the speed of sound"
        " in the plasma,\"\n"
        "    Nucl. Phys. B 891 (2015), 159-199 [arXiv:1410.3875 [hep-ph]].\n"
        "  [Espinosa+ (2010)]\n"
        "    J. R. Espinosa and J. M. No,\n"
        "    \"Energy Budget of Cosmological First-order Phase Transitions,\"\n"
        "    JCAP 06 (2010), 028 [arXiv:1004.4187 [hep-ph]].";

    m.attr("__version__") = HYDRODEX_VERSION_STRING;

    // replace GSL's default error handler (abort, which kills the interpreter)
    // so that GSL errors are raised as exceptions, but keep any error handler
    // that another library has already installed
    gsl_error_handler_t* previous_handler = gsl_set_error_handler_off();
    if ( previous_handler != nullptr ) {
        gsl_set_error_handler(previous_handler);
    }

    bind_settings(m);
    bind_utils(m);

    py::module_ pcalc = m.def_submodule(
        "profile_calculator", 
        "Functions to calculate the fluid profile of a bubble in a cosmological"
        " phase transition."
    );
    bind_profile_calculator(pcalc);
    py::module_ ssmodel = m.def_submodule(
        "sound_shell_model", 
        "Functions to calculate the gravitational wave spectrum within the"
        " sound shell model."
    );
    bind_sound_shell_model(ssmodel);

    bind_fluid_profile(m);
    bind_sound_shell_spectrum(m);

} // PYBIND11_MODULE(hydrodex, m)
