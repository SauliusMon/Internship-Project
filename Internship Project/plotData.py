import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

fig, ax = plt.subplots()

# Use a raw string (r"...") to prevent Windows backslashes from acting as escape characters
file_path = r"C:\Users\Sauls\source\repos\Internship Project\Internship Project\tof_rezultatai.csv"

# 1. Initialize the plot layout and an empty line object once
line, = ax.plot([], [], color='blue', marker='o', linewidth=0.5, markersize=3)
ax.set_title("TOF")
ax.set_xlabel("t")
ax.set_ylabel("I")

def atnaujinti_grafika(i):
    try:
        # 2. Check if file exists and isn't empty before trying to load it
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            return

        # 3. Read directly with numpy (avoids the slow double-read of f.readlines())
        data = np.loadtxt(file_path, delimiter=',', skiprows=1)
        
        if data.size == 0:
            return
            
        if data.ndim == 1:
            t = np.array([data[0]])
            I_total = np.array([data[1]])
        else:
            t = data[:, 0]
            I_total = data[:, 1]

        # 4. Update the existing line data instead of clearing the whole axis
        line.set_data(t, I_total)
        
        # 5. Rescale the axes to fit the new incoming data dynamically
        ax.relim()
        ax.autoscale_view()

        if len(t) > 1:
            charge_Q = np.trapz(I_total, x=t) 
            ax.set_title(f"TOF | Integral (Charge Q): {charge_Q:.4e}")
        
    except Exception as e:
        # Don't silently pass - if the file is locked by the writing program, 
        # it's better to know why it skipped a frame
        pass # Replaced with pass to keep your original behavior, but print(e) is recommended for debugging

ani = animation.FuncAnimation(fig, atnaujinti_grafika, interval=100, cache_frame_data=False)
plt.show()