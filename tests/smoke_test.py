"""Smoke test of the HydróDex Python interface.

Checks that the compiled module can be imported and that a fluid profile and a
(coarsely sampled) gravitational wave spectrum can be computed and give finite
results. This is not a test of the numerical accuracy.

Usage:
    python tests/smoke_test.py [BUILD_DIR]

The module is loaded from BUILD_DIR (default: ./build), so that the freshly
built module takes precedence over any installed version.
"""

import glob
import importlib.util
import math
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_module(build_dir):
    pattern = os.path.join(
        build_dir, "lib", "python*", "site-packages", "hydrodex*.so"
    )
    candidates = glob.glob(pattern)
    if len(candidates) != 1:
        sys.exit(f"Expected exactly one module matching {pattern}, "
                 f"found {candidates}")
    spec = importlib.util.spec_from_file_location("hydrodex", candidates[0])
    module = importlib.util.module_from_spec(spec)
    sys.modules["hydrodex"] = module
    spec.loader.exec_module(module)
    return module


def project_version():
    with open(os.path.join(ROOT, "CMakeLists.txt")) as f:
        match = re.search(r"project\([^)]*VERSION\s+([0-9.]+)", f.read())
    return match.group(1)


def check(condition, message):
    if not condition:
        sys.exit(f"FAILED: {message}")
    print(f"ok: {message}")


def all_finite(values):
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def main():
    build_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build")
    hydrodex = load_module(build_dir)
    pc = hydrodex.profile_calculator
    ssm = hydrodex.sound_shell_model

    check(hydrodex.__version__ == project_version(),
          f"__version__ ({hydrodex.__version__}) matches CMakeLists.txt")

    # fluid profile
    profile = hydrodex.FluidProfile(xi_w=0.65, alpha_N=0.05)
    xi = np.asarray(profile.get_xi_coordinate(), dtype=float)
    v = np.asarray(profile.get_fluid_velocity(), dtype=float)
    check(xi.size > 0 and xi.size == v.size, "profile arrays are non-empty")
    check(all_finite(xi) and all_finite(v), "profile arrays are finite")
    check(np.all((v >= 0.0) & (v < 1.0)), "fluid velocity in [0, 1)")
    kappa = profile.get_efficiency(pc.EfficiencyNormalization.TotalEnergyDensity)
    check(math.isfinite(kappa) and 0.0 < kappa < 1.0,
          f"efficiency in (0, 1): {kappa}")

    # gravitational wave spectrum (coarse grids, results not converged)
    spec = hydrodex.SoundShellSpectrum(
        profile=profile,
        H_R=1.0e-2,
        nhist=ssm.NucleationHistory.EXPONENTIAL,
        exhist=ssm.ExpansionHistory.RADIATION_DOMINATION,
    )
    spec.set_source_duration(1.0e-1)
    spec.set_correlator_integration_settings(20, 16)
    spec.calculate_A_squared(1.0e-2, 1.0e3, 500)
    spec.calculate_kinetic_spectrum(1.0e-2, 1.0e3, 100)
    spec.calculate_gw_spectrum(1.0e-1, 1.0e2, 10, 100)

    Uf2 = spec.get_mean_squared_fluid_velocity()
    check(math.isfinite(Uf2) and Uf2 > 0.0,
          f"mean squared fluid velocity positive: {Uf2}")
    Omega_GW = np.asarray(spec.get_gravitational_wave_spectrum(), dtype=float)
    check(Omega_GW.size == 10, "GW spectrum has the requested size")
    check(all_finite(Omega_GW) and np.all(Omega_GW >= 0.0),
          "GW spectrum is finite and non-negative")

    print("All smoke tests passed.")


if __name__ == "__main__":
    main()
