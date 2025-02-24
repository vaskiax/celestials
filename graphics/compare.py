import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from matplotlib import rcParams

def animations(r1, r2, r3, save=False, fps=1.10):
    """
    Smooth animation of the system motion with optimized performance.
    Fixes issues where rerunning the code accumulates choppiness.
    """

    # 🔹 Reset Matplotlib to avoid frame accumulation
    plt.close('all')  

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter3D(r1[:, 0], r1[:, 1], r1[:, 2], 
                 c='black', marker='o', s=0.05, alpha=0.5)

    point1, = ax.plot([], [], [], 'ro', markersize=8, label='Guided Center')  
    point2, = ax.plot([], [], [], 'bo', markersize=8, label='Keplerian')  

    ax.set_axis_off()
    ax.view_init(elev=90, azim=0)
    plt.legend(loc='upper right')
    plt.title('Comparison: Guided Center vs Keplerian Orbits')

    def update(frame):
        """Update function for animation"""
        x1, y1, z1 = r2[frame, 0], r2[frame, 1], r2[frame, 2]
        x2, y2, z2 = r3[frame, 0], r3[frame, 1], r3[frame, 2]

        point1.set_data([x1], [y1])
        point1.set_3d_properties([z1])

        point2.set_data([x2], [y2])
        point2.set_3d_properties([z2])

        return point1, point2  

    # 🔹 Prevent frame accumulation by setting a fixed frame count
    frames = min(len(r1), len(r2), len(r3), 720)
    interval = 100/fps  # 30 FPS
    rcParams['animation.embed_limit'] = 2**128

    ani = FuncAnimation(fig, update, frames=frames, interval=interval, cache_frame_data=False)

    if save:
        plt.close(fig) 
        ani.save('media/comparison.gif', writer='ffmpeg', fps=60)

    plt.close(fig) 
    return HTML(ani.to_jshtml())  
