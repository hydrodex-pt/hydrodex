# Changelog

All notable changes to HydróDex are documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.0] - 2026-11-XX

Initial release of HydróDex, accompanying the preprint "HydróDex: gotta match
'em all!" by G. Barni and E. Madge.

HydróDex is the successor of
[inverse_pt](https://github.com/eric-madge/inverse_pt). Major changes compared
to inverse_pt:
- Fluid profiles in the constant-sound-speed model with different speeds of
  sound in the two phases, instead of only the bag model
- Gravitational wave spectrum for a universe with a constant
  equation-of-state parameter, in addition to radiation domination and a flat
  (non-expanding) universe

### Added
- C++ library `libhydrodex` for computing the fluid profiles of direct and
  inverse cosmological phase transitions (`hydrodex::profile_calculator`,
  `hydrodex::FluidProfile`)
- Computation of the gravitational wave spectrum within the sound shell model
  (`hydrodex::sound_shell_model`, `hydrodex::SoundShellSpectrum`)
- Python interface `hydrodex`, installable via CMake or pip
- Command-line executable `hydrodex`
- Doxygen documentation of the C++ library and Sphinx documentation of the
  Python interface
- Example Jupyter notebooks, including the code to reproduce the figures of
  the accompanying paper
