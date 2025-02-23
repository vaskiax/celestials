def animations(r1, r2, r3, save=False):
    """
    Smooth animation of the system motion with optimized performance.
    """
    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation
    from IPython.display import HTML
    from matplotlib import rcParams

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Background scatter points
    ax.scatter3D(r1[:, 0], r1[:, 1], r1[:, 2], 
                 c='black', marker='o', s=0.05, alpha=0.5)

    # Create TWO separate animated points
    point1, = ax.plot([], [], [], 'ro', markersize=8, label='Guided Center')  # Red point
    point2, = ax.plot([], [], [], 'bo', markersize=8, label='Keplerian')  # Blue point
    ax.set_axis_off()
    ax.view_init(elev=90, azim=0)
    plt.legend()
    plt.title('Keplerian vs Guided Center approach comparison')
    

    def update(frame):
        """Update function for animation"""
        x1, y1, z1 = r2[frame, 0], r2[frame, 1], r2[frame, 2]
        x2, y2, z2 = r3[frame, 0], r3[frame, 1], r3[frame, 2]

        # Update first point
        point1.set_data([x1], [y1])
        point1.set_3d_properties([z1])

        # Update second point
        point2.set_data([x2], [y2])
        point2.set_3d_properties([z2])

        return point1, point2  # Must return both!

    frames = len(r1)
    interval = 1000 / 30  # 🔥 30 FPS → 33ms per frame
    rcParams['animation.embed_limit'] = 2**128
    ani = FuncAnimation(fig, update, frames=frames, interval=interval, blit=True)

    if save:
        ani.save('media/comparison.gif', writer='imagemagick', fps=30)

    plt.close()
    return HTML(ani.to_jshtml())
