import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np



U = 3.0

eps_0 = 8.854e-12;       
eps = 3.0;                
d = 1e-6;                  
S = 4e-6;                  

C =(eps_0 * eps * S) / d

# Sukuriame 2x2 grafikų langą
fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, figsize=(12, 8))

base_path = r"C:\Users\Sauls\source\repos\Internship Project\Internship Project\\"
file_results = base_path + "tof_rezultatai.csv"
file_profile = base_path + "tof_E_profilis.csv"

# --- 1. Kairė Viršuje: Srovė I(t) ---
(line_I,) = ax1.plot([], [], color="blue", linewidth=1.2)
ax1.set_title("TOF - Srovė I(t)")
ax1.set_ylabel("I (A)")
ax1.set_xlabel("t (s)")
ax1.grid(True)

# --- 2. Kairė Apačioje: Integral E(t) ---
(line_E_int,) = ax2.plot([], [], color="red", linewidth=1.2)
ax2.set_title("Lauko integralas laike E_int(t)")
ax2.set_ylabel("Integral E (V)")
ax2.set_xlabel("t (s)")
ax2.grid(True)

# --- 3. Dešinė Viršuje: E laukas nuo koordinatės E(x) ---
(line_E_prof,) = ax3.plot([], [], color="green", linewidth=1.5)
ax3.set_title("E laukas nuo koordinatės E(x)")
ax3.set_ylabel("E (V/m)")
ax3.set_xlabel("x (m)")
ax3.grid(True)

# --- 4. Dešinė Apačioje: n(x) ir p(x) koncentracijos ---
(line_n,) = ax4.plot([], [], color="cyan", linewidth=1.5, label="n (elektronai)")
(line_p,) = ax4.plot([], [], color="magenta", linewidth=1.5, label="p (skylės)")
ax4.set_title("Krūvininkų koncentracija n(x) ir p(x)")
ax4.set_ylabel("Koncentracija (1/m³)")
ax4.set_xlabel("x (m)")
ax4.legend(loc="upper right")
ax4.grid(True)


ax3.set_xlim(0, d);
ax4.set_xlim(0, d);


def atnaujinti_grafika(i):
    try:
        # --- Laiko duomenys (ax1, ax2) ---
        if os.path.exists(file_results) and os.path.getsize(file_results) > 0:
            data = np.loadtxt(file_results, delimiter=",", skiprows=1)
            if data.size > 0:
                if data.ndim == 1:
                    t, I_total, integral_E = data[0], data[1], data[2]
                    t, I_total, integral_E = (
                        np.array([t]),
                        np.array([I_total]),
                        np.array([integral_E]),
                    )
                else:
                    t, I_total, integral_E = data[:, 0], data[:, 1], data[:, 2]

                line_I.set_data(t, I_total)
                line_E_int.set_data(t, integral_E)

                ax1.relim()
                ax1.autoscale_view()
                ax2.relim()
                ax2.autoscale_view()

                if len(t) > 1:
                    charge_Q = (
                        np.trapezoid(I_total, x=t)
                        if hasattr(np, "trapezoid")
                        else np.trapz(I_total, x=t)
                    )
                    L = charge_Q / C / U
                    
                    ax1.set_title(
                        f"TOF | Srovė I(t) | L: {L:.4e}"
                    )

        # --- Spatialiniai duomenys (ax3, ax4) ---
        if os.path.exists(file_profile) and os.path.getsize(file_profile) > 0:
            prof_data = np.loadtxt(file_profile, delimiter=",", skiprows=1)
            if prof_data.size > 0 and prof_data.ndim == 2:
                x_coords = prof_data[:, 0]
                E_values = prof_data[:, 1]
                n_values = prof_data[:, 2]
                p_values = prof_data[:, 3]

                # Atnaujiname E(x)
                line_E_prof.set_data(x_coords, E_values)
                ax3.relim()
                ax3.autoscale_view()

                # Atnaujiname n(x) ir p(x)
                line_n.set_data(x_coords, n_values)
                line_p.set_data(x_coords, p_values)
                ax4.relim()
                ax4.autoscale_view()

    except Exception:
        pass


ani = animation.FuncAnimation(
    fig, atnaujinti_grafika, interval=100, cache_frame_data=False
)

plt.tight_layout()
plt.show()