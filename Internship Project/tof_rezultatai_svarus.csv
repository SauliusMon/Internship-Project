import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

fig, ax = plt.subplots()
file_path = "C:\\Users\\Sauls\\source\\repos\\Internship Project\\Internship Project\\tof_rezultatai.csv"

def atnaujinti_grafika(i):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        if len(lines) <= 1:
            return

        data = np.loadtxt(file_path, delimiter=',', skiprows=1)
        
        if data.ndim == 1 and data.size > 0:
            t = np.array([data[0]])
            I_total = np.array([data[1]])
        elif data.ndim == 2 and data.shape[0] > 0:
            t = data[:, 0]
            I_total = data[:, 1]
        else:
            return

        ax.clear()
        ax.plot(t, I_total, color='blue', marker='o', linewidth=0.5, markersize=3)
        ax.set_title("TOF")
        ax.set_xlabel("t")
        ax.set_ylabel("I")
        
    except Exception:
        pass

ani = animation.FuncAnimation(fig, atnaujinti_grafika, interval=100, cache_frame_data=False)
plt.show()