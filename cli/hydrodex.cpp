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

#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

#include "hydrodex/fluid_profile.hpp"
#include "hydrodex/sound_shell_spectrum.hpp"
#include "hydrodex/profile_calculator.hpp"
#include "hydrodex/sound_shell_model.hpp"
#include "hydrodex/settings.hpp"
#include "hydrodex/constants.hpp"
#include "hydrodex/utils.hpp"

using namespace hydrodex::utils;
using namespace hydrodex::profile_calculator;
using namespace hydrodex::sound_shell_model;
using namespace hydrodex::constants;

int main(int argc, char* argv[]) {

    hydrodex::settings::set_gsl_abort_on_error(false);
    //hydrodex::settings::debug = true;
    //hydrodex::settings::consistency_check=hydrodex::settings::CheckMode::ERROR;

    if (argc < 3) {
        std::cerr << "Usage: ./hydrodex <vw> <alpha_N> [HR] [Hdtau] [c_s]"
                  << " [omega] [Omega_star] [c_plus] [c_minus] [precision]\n"
                  << "Example: ./hydrodex"
                  << " 0.4 0.1 0.01 1.0 0.577 0.667 1.0 0.577 0.577 0.001\n";
        return 1;
    }
    double vw = std::stod(argv[1]);
    double alpha_N = std::stod(argv[2]);
    double H_R = (argc >= 4) ? std::stod(argv[3]) : 1.0;
    double H_delta_tau = (argc >= 5) ? std::stod(argv[4]) : 1.0;
    double c_s = (argc >= 6) ? std::stod(argv[5]) : cs_bag;
    double omega = (argc >= 7) ? std::stod(argv[6]) : 1.0/3.0;
    double Omega_star = (argc >= 8) ? std::stod(argv[7]) : 1.0;
    double c_plus = (argc >= 9) ? std::stod(argv[8]) : cs_bag;
    double c_minus = (argc >= 10) ? std::stod(argv[9]) : cs_bag;
    double precision = (argc >= 11) ? std::stod(argv[10]) : 1e-3;

    std::chrono::time_point<std::chrono::high_resolution_clock> main_start 
        = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> sub_start;
    std::chrono::duration<double> duration;

    std::vector<double> xi, v;

    std::cout << "Calculating fluid profile for xi_w = " << vw 
              << " and alpha_N = " << alpha_N << " (c+ = " << c_plus
              << ", c- = " << c_minus << ")." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();

    HydroSolutionType sol_type = classify_wave_fluid_frame(
        vw, alpha_N, c_plus, c_minus, precision
    );
    double alpha_plus = find_alpha_plus(
        sol_type, alpha_N, vw, c_plus, c_minus, precision
    );
    auto [vplus, vminus] = compute_fluid_velocities(
        vw, alpha_plus, sol_type, c_plus, c_minus
    );

    std::cout << "Found " << sol_type << " solution with alpha_+ = " 
              << alpha_plus << ".\nFluid velocities: v_+ = " << vplus
              << ", v_- = " << vminus << ".\nCalculating velocity profile."
              << std::endl;

    size_t idx_wall = compute_velocity_profile(
        sol_type, xi, v, vw, alpha_plus, c_plus, c_minus, precision
    );
    size_t idx_shock = NO_INDEX;
    if ( 
        sol_type == HydroSolutionType::Deflagration 
        || sol_type == HydroSolutionType::Hybrid 
        || sol_type == HydroSolutionType::HybridPlus
        || sol_type == HydroSolutionType::InverseDeflagration
        || sol_type == HydroSolutionType::InverseHybrid
    ) {
        idx_shock = get_shock_index(sol_type, xi.size());
    }
    std::cout << "Velocity profile calculated with " << xi.size() << " points."
              << "\n  Velocity jumps from " << v[idx_wall-1] << " to " 
              << v[idx_wall] << " at the bubble wall, xi = " << xi[idx_wall] 
              << " (index: " << idx_wall << ")." ;
    if ( idx_shock != NO_INDEX ) { 
        std::cout << "\n  Velocity jumps from " << v[idx_shock-1] << " to " 
                  << v[idx_shock] << " at the shock position, xi = " 
                  << xi[idx_shock] << " (index: " << idx_shock << ").";
    }
    std::cout << std::endl;
    std::vector<double> w = compute_enthalpy_density(
        xi, v, idx_wall, idx_shock, c_plus, c_minus
    );
    double wN = w.back();
    std::cout << "Enthalpy density calculated. wN: " << wN << std::endl;

    double epsilon = compute_epsilon(wN, alpha_N, c_plus, c_minus);
    std::cout << "Vacuum energy: epsilon = " << epsilon << std::endl;
    std::vector<double> e = compute_energy_density_profile(
        w, idx_wall, epsilon, c_plus, c_minus
    );
    std::vector<double> p = compute_pressure_profile(
        w, idx_wall, epsilon, c_plus, c_minus
    );

    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "All profiles calculated after " << duration.count() 
              << " seconds. Writing to file." << std::endl;
    std::ofstream outfile("profiles.csv");
    outfile << "# xi, v, w, e, p\n";
    for ( size_t i=0; i<xi.size(); ++i ) {
        outfile << xi[i] << ", " << v[i] << ", " << w[i] << ", " << e[i] << ", "
                << p[i] << "\n";
    }
    outfile.close();

    hydrodex::FluidProfile::set_verbose(true);
    hydrodex::FluidProfile fluid_profile(
        xi, v, w, e, p, idx_wall, idx_shock, c_plus, c_minus
    );

    std::cout << "Calculating GW spectrum with H*R = " << H_R 
              << " and H*deltaTau = " << H_delta_tau << ":" << std::endl;

    hydrodex::SoundShellSpectrum::set_verbose(true);
    hydrodex::SoundShellSpectrum spectra(
        fluid_profile, 
        H_R,
        NucleationHistory::EXPONENTIAL, 
        ExpansionHistory::RADIATION_DOMINATION,
        c_s,
        omega,
        Omega_star
    );

    spectra.set_chi(1e-2, 1e4, 2401); 
    spectra.set_kR(1e-3, 1e4, 200);
    spectra.set_qR_from_kR(10001);

    std::vector<double> chi = spectra.get_chi();
    std::vector<double> qR = spectra.get_qR();
    std::vector<double> kR = spectra.get_kR();

    if ( H_delta_tau > 0. ) { spectra.set_source_duration(H_delta_tau); }
    size_t N_p = 200, N_z = 64;
    spectra.set_correlator_integration_settings(N_p, N_z);

    sub_start = std::chrono::high_resolution_clock::now();
    auto [df, l] = spectra.calculate_A_squared();
    std::vector<double> A2 = spectra.get_shape_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Shape function calculated after " << duration.count() 
              << " seconds." << std::endl;

    spectra.set_nucleation_history(NucleationHistory::EXPONENTIAL);
    std::cout << "Calculating spectra for exponential nucleation:" << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> E_kin_over_R_exp = spectra.get_kinetic_spectrum();
    double K_peak_kin_exp = spectra.get_kinetic_spectrum_peak_position();
    double E_kin_star_over_R_exp = spectra.get_kinetic_spectrum_amplitude();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic spectrum calculated after " << duration.count() 
              << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double Uf2_exp = spectra.get_mean_squared_fluid_velocity();
    double Omega_K_exp = spectra.get_kinetic_energy_density_parameter();
    double curly_K_exp = spectra.get_normalized_kinetic_density_parameter();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic density parameter calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> k_E_Pi_exp = 
        spectra.get_anisotropic_stress_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Anisotropic stress power spectrum calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double curly_C_exp = spectra.get_anisotropic_stress_UETC_normalization();
    std::vector<double> K3_zeta_Pi_exp = 
        spectra.get_normalized_anisotropic_stress_power_spectrum();
    auto [K_peak_Pi_exp, K3_zeta_Pi_peak_exp] = find_peak(
        spectra.get_kR(), K3_zeta_Pi_exp
    );
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized anisotropic stress calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::FLAT_UNIVERSE);
    double tilde_Delta_0_flat_exp = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::INFINITE_DURATION
    );
    std::vector<double> Omega_GW_approx_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated approximate GW spectrum after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::FULL
    );
    std::vector<double> Omega_GW_flat_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_flat_exp = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_flat_exp =
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::RADIATION_DOMINATION);
    double tilde_Delta_0_rad_exp = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a radiation dominated"
              << " universe after " << duration.count() << " seconds."
              << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> Omega_GW_rad_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a radiation dominated universe"
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_rad_exp = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_rad_exp = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a radiation dominated universe"
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::CONSTANT_EOS);
    double tilde_Delta_0_omega_exp = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a universe with omega"
              << " = " << omega << " after " << duration.count() << " seconds."
              << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> Omega_GW_omega_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a universe with omega = " << omega 
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_omega_exp = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_omega_exp = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a universe with omega = " << omega
              << " after " << duration.count() << " seconds." << std::endl;
    spectra.set_expansion_history(ExpansionHistory::RADIATION_DOMINATION);

    spectra.set_nucleation_history(NucleationHistory::SIMULTANEOUS);
    std::cout << "Calculating spectra for simultaneous nucleation:" 
              << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> E_kin_over_R_sim = spectra.get_kinetic_spectrum();
    double K_peak_kin_sim = spectra.get_kinetic_spectrum_peak_position();
    double E_kin_star_over_R_sim = spectra.get_kinetic_spectrum_amplitude();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic spectrum calculated after " << duration.count() 
              << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double Uf2_sim = spectra.get_mean_squared_fluid_velocity();
    double Omega_K_sim = spectra.get_kinetic_energy_density_parameter();
    double curly_K_sim = spectra.get_normalized_kinetic_density_parameter();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic density parameter calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> k_E_Pi_sim = 
        spectra.get_anisotropic_stress_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Anisotropic stress power spectrum calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double curly_C_sim = spectra.get_anisotropic_stress_UETC_normalization();
    std::vector<double> K3_zeta_Pi_sim =
        spectra.get_normalized_anisotropic_stress_power_spectrum();
    auto [K_peak_Pi_sim, K3_zeta_Pi_peak_sim] = find_peak(
        spectra.get_kR(), K3_zeta_Pi_sim
    );
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized anisotropic stress calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::FLAT_UNIVERSE);
    double tilde_Delta_0_flat_sim = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a flat universe after "
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::INFINITE_DURATION
    );
    std::vector<double> Omega_GW_approx_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated approximate GW spectrum after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::FULL
    );
    std::vector<double> Omega_GW_flat_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_flat_sim = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_flat_sim = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::RADIATION_DOMINATION);
    double tilde_Delta_0_rad_sim = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a radiation dominated"
              << " universe after " << duration.count() << " seconds."
              << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> Omega_GW_rad_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a radiation dominated universe"
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_rad_sim = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_rad_sim = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a radiation dominated universe"
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::CONSTANT_EOS);
    double tilde_Delta_0_omega_sim = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a universe with omega"
              << " = " << omega << " after " << duration.count() << " seconds."
              << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> Omega_GW_omega_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a universe with omega = " << omega 
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_omega_sim = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_omega_sim = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a universe with omega = " << omega
              << " after " << duration.count() << " seconds." << std::endl;
    spectra.set_expansion_history(ExpansionHistory::RADIATION_DOMINATION);

    duration = std::chrono::high_resolution_clock::now() - main_start;
    std::cout << "\nTotal run time: " << duration.count() << " seconds.\n"
              << "Results:\n"
              << "  exponential nucleation:\n"
              << "    10^4 E_kin_*/R_* = " << 1.0e4*E_kin_star_over_R_exp <<"\n"
              << "    K_kin_peak       = " << K_peak_kin_exp << "\n"
              << "    curly_K          = " << curly_K_exp << "\n"
              << "    200 Uf^2         = " << 2.0e2*Uf2_exp << "\n"
              << "    10^2 Omega_K     = " << 1.0e2*Omega_K_exp << "\n"
              << "    curly_C          = " << curly_C_exp << "\n"
              << "    K_Pi_peak        = " << K_peak_Pi_exp << "\n"
              << "    K^3 zeta_Pi_peak = " << K3_zeta_Pi_peak_exp << "\n"
              << "    Delta~0 (rad.)   = " << tilde_Delta_0_rad_exp << "\n"
              << "    Delta~0 (stat.)  = " << tilde_Delta_0_flat_exp << "\n"
              << "    Delta~0 (omega)  = " << tilde_Delta_0_omega_exp << "\n"
              << "  simultaneous nucleation:\n"
              << "    10^4 E_kin_*/R_* = " << 1.0e4*E_kin_star_over_R_sim <<"\n"
              << "    K_kin_peak       = " << K_peak_kin_sim << "\n"
              << "    curly_K          = " << curly_K_sim << "\n"
              << "    200 Uf^2         = " << 2.0e2*Uf2_sim << "\n"
              << "    10^2 Omega_K     = " << 1.0e2*Omega_K_sim << "\n"
              << "    curly_C          = " << curly_C_sim << "\n"
              << "    K_Pi_peak        = " << K_peak_Pi_sim << "\n"
              << "    K^3 zeta_Pi_peak = " << K3_zeta_Pi_peak_sim << "\n"
              << "    Delta~0 (rad.)   = " << tilde_Delta_0_rad_sim << "\n"
              << "    Delta~0 (stat.)  = " << tilde_Delta_0_flat_sim << "\n"
              << "    Delta~0 (omega)  = " << tilde_Delta_0_omega_sim << "\n"
              << std::endl;

    std::cout << "Done. Writing spectra to file" << std::endl;
    outfile.open("spectra.csv");
    outfile << "# chi, df_dchi, l, A2, qR, Ekin_exp, Ekin_sim, kR, kEPi_exp, "
            << "kEPi_sim, OmegaGW_approx_exp, OmegaGW_approx_sim, "
            << "OmegaGW_flat_exp, OmegaGW_flat_sim, OmegaGW_exp, OmegaGW_sim,"
            << " OmegaGW_omega_exp, OmegaGW_omega_sim\n";
    for ( size_t i=0; i<kR.size() || i<qR.size() || i<chi.size(); ++i ) {
        if ( i < chi.size() ) { 
            outfile << chi[i] << ", " << df[i] << ", " << l[i] << ", "
                    << A2[i] << ", "; 
        }
        else { outfile << ",,,, "; }
        if ( i < qR.size() ) { 
            outfile << qR[i] << ", " << E_kin_over_R_exp[i] << ", " 
                    << E_kin_over_R_sim[i] << ", "; 
        }
        else { outfile << ",,, "; }
        if ( i < kR.size() ) { 
            outfile << kR[i] << ", " << k_E_Pi_exp[i]
                             << ", " << k_E_Pi_sim[i] 
                             << ", " << Omega_GW_approx_exp[i] 
                             << ", " << Omega_GW_approx_sim[i] 
                             << ", " << Omega_GW_flat_exp[i] 
                             << ", " << Omega_GW_flat_sim[i] 
                             << ", " << Omega_GW_rad_exp[i] 
                             << ", " << Omega_GW_rad_sim[i] 
                             << ", " << Omega_GW_omega_exp[i] 
                             << ", " << Omega_GW_omega_sim[i] 
                             << "\n"; 
        } else { outfile << ",,,,,,,,,,\n"; }
    }
    outfile.close();
    
    std::cout << "Done." << std::endl;

    return 0;
}