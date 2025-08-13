import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d
from libra_py.units import au2ev, au2fs
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Bi-exponential decay function with fixed E_inf
# ----------------------------
def bi_exponential_decay(t, A1, tau1, A2, tau2, Einf):
    return Einf + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)

# ----------------------------
# Load Adiabatic Energies
# ----------------------------
def load_adiabatic_energies(start_step, end_step, path_template, scale=au2ev, time_scale=au2fs):
    energies = []
    time = None
    for step in range(start_step, end_step):
        path = path_template.format(step=step)
        energy = sp.load_npz(path).todense().real
        if time is None:
            time = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time

# ----------------------------
# Main Analysis Function
# ----------------------------
def calculate_total_energy_decay_to_S1(
    adiabatic_energies, methods, icond_range, method_config,
    smoothing_window=100, fit_time_cutoff=1000, verbose=True
):
    energy_results = {}
    average_energies = {}
    time_results = {}
    decay_fit_results = {}

    all_pops = []
    all_E1 = []

    # --- Compute weighted average S1 energy ---
    for method in methods:
        path_template = method_config[method]["path_template"]
        for icond in icond_range:
            with h5py.File(path_template.format(icond=icond), 'r') as f:
                sh_pop = np.array(f['sh_pop_adi/data'])
                n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                pops = sh_pop[:n_steps, 1]
                E1 = adiabatic_energies[:n_steps, 1]
                all_pops.append(pops)
                all_E1.append(E1)

    all_pops = np.array(all_pops)
    all_E1 = np.array(all_E1)
    E1_avg = np.sum(all_pops * all_E1) / np.sum(all_pops)

    for method in methods:
        print(f"\nProcessing method: {method}")
        try:
            method_results = []
            method_times = None
            path_template = method_config[method]["path_template"]

            tau1_list = []
            tau2_list = []
            initial_total_energies = []

            for icond in icond_range:
                path = path_template.format(icond=icond)
                with h5py.File(path, 'r') as f:
                    sh_pop = np.array(f['sh_pop_adi/data'])
                    time = np.array(f['time/data']) * au2fs

                    n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                    sh_pop = sh_pop[:n_steps, :]
                    time = time[:n_steps]
                    Hvib_rolled = np.roll(adiabatic_energies[:n_steps, :], -icond, axis=0)

                    total_energy = np.sum(sh_pop * Hvib_rolled, axis=1)
                    initial_total_energies.append(total_energy[0])

                    # Smooth individual energy
                    smoothed = uniform_filter1d(total_energy, size=smoothing_window)

                    # Fit to bi-exponential for this trajectory
                    fit_mask = time < fit_time_cutoff
                    fit_times = time[fit_mask]
                    fit_values = smoothed[fit_mask]

                    def fit_func_biexp(t, A1, tau1, A2, tau2):
                        return bi_exponential_decay(t, A1, tau1, A2, tau2, E1_avg)

                    E0_i = total_energy[0]
                    p0 = [E0_i - E1_avg, 300, (E0_i - E1_avg) * 0.5, 1000]
                    bounds = ([-np.inf, 1e-2, -np.inf, 1e-2], [np.inf, 1e4, np.inf, 1e4])

                    try:
                        popt, _ = curve_fit(fit_func_biexp, fit_times, fit_values, p0=p0, bounds=bounds)
                        A1_i, tau1_i, A2_i, tau2_i = popt
                        tau1_list.append(tau1_i)
                        tau2_list.append(tau2_i)
                    except Exception as e_fit:
                        print(f"  Fit failed for icond {icond}: {e_fit}")

                    method_results.append(total_energy)
                    if method_times is None:
                        method_times = time

            # Aggregate energy and fit summary
            total_energy_avg = np.mean(method_results, axis=0)
            total_energy_std = np.std(method_results, axis=0)
            smoothed_total = uniform_filter1d(total_energy_avg, size=smoothing_window)

            avg_tau1 = np.mean(tau1_list)
            std_tau1 = np.std(tau1_list)
            avg_tau2 = np.mean(tau2_list)
            std_tau2 = np.std(tau2_list)

            # Fit average energy for plotting (dashed black line)
            fit_mask = method_times < fit_time_cutoff
            fit_times = method_times[fit_mask]
            fit_values = smoothed_total[fit_mask]

            def fit_func_biexp(t, A1, tau1, A2, tau2):
                return bi_exponential_decay(t, A1, tau1, A2, tau2, E1_avg)

            E0_fixed = np.mean(initial_total_energies)
            p0 = [E0_fixed - E1_avg, 300, (E0_fixed - E1_avg) * 0.5, 1000]

            try:
                popt, _ = curve_fit(fit_func_biexp, fit_times, fit_values, p0=p0,
                                    bounds=([-np.inf, 1e-2, -np.inf, 1e-2], [np.inf, 1e4, np.inf, 1e4]))
                fit_curve = fit_func_biexp(method_times, *popt)
                tau1_fit, tau2_fit = popt[1], popt[3]
            except Exception as e_avg_fit:
                print(f"Average fit failed: {e_avg_fit}")
                fit_curve = None

            # Plot
            if verbose:
                plt.figure()
                plt.plot(method_times, total_energy_avg, label='Raw Total Energy', alpha=0.4)
                plt.fill_between(
                    method_times,
                    total_energy_avg - total_energy_std,
                    total_energy_avg + total_energy_std,
                    color='gray', alpha=0.3, label='±1σ'
                )
                plt.plot(method_times, smoothed_total, label='Smoothed Avg Energy', color='red')
                if fit_curve is not None:
                    plt.plot(method_times, fit_curve, '--', color='black',
                             label=f'Avg Fit: τ₁={tau1_fit:.1f} fs, τ₂={tau2_fit:.1f} fs')

                plt.axhline(E0_fixed, color='green', linestyle=':', label=f'Initial ≈ {E0_fixed:.3f} eV')
                plt.axhline(E1_avg, color='blue', linestyle=':', label=f'S₁ Energy ≈ {E1_avg:.3f} eV')

                plt.xlabel("Time (fs)", fontsize=22)
                plt.ylabel("Average Excitation Energy (eV)", fontsize=22)
                plt.tick_params(axis='both', which='major', labelsize=22)
                plt.gcf().set_size_inches(12, 6)
                plt.tight_layout(rect=[0, 0, 0.75, 1])
              #  plt.legend(fontsize=14)
                plt.show()

            print(f"τ₁ = {avg_tau1:.2f} ± {std_tau1:.2f} fs")
            print(f"τ₂ = {avg_tau2:.2f} ± {std_tau2:.2f} fs")

            energy_results[method] = (method_results, (avg_tau1, std_tau1, avg_tau2, std_tau2))
            average_energies[method] = smoothed_total
            time_results[method] = method_times
            decay_fit_results[method] = {
                "tau1_avg": avg_tau1,
                "tau1_std": std_tau1,
                "tau2_avg": avg_tau2,
                "tau2_std": std_tau2,
                "Einf": E1_avg,
                "individual_tau1": tau1_list,
                "individual_tau2": tau2_list
            }

        except Exception as e:
            print(f"Error in method {method}: {e}")

    return energy_results, average_energies, time_results, decay_fit_results


# ----------------------------
# Energy file info
# ----------------------------
start_step = 1000
end_step = 3990
energy_path_template = "../step3_NACs/results/Hvib_ci_{step}_re.npz"
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Load adiabatic energies
# ----------------------------
adiabatic_energies, _ = load_adiabatic_energies(start_step, end_step, energy_path_template)

# ----------------------------
# Method configuration
# ----------------------------
method_config = {
    "DISH": {
        "path_template": "/DISH_results/DISH_icond_{icond}/mem_data.hdf"
    }
}

# ----------------------------
# Initial condition indices
# ----------------------------
icond_range = range(1, 3000, 100)

# ----------------------------
# Run energy decay analysis with fixed decay to S₁
# ----------------------------
methods = ["DISH"]
total_results, total_avg, total_times, total_fits = calculate_total_energy_decay_to_S1(
    adiabatic_energies=adiabatic_energies,
    methods=methods,
    icond_range=icond_range,
    method_config=method_config,
    smoothing_window=100,
    verbose=True
)