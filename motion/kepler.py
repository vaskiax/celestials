class Kepler:
    """
    Kepler's equation solver.

    This class provides a method to solve Kepler's equation for the eccentric
    anomaly E given the mean anomaly M and the eccentricity e.

    Attributes
    ----------
    M : float
        The mean anomaly in radians.
    e : float
        The orbital eccentricity.
    E : float    
        The eccentric anomaly in radians.
    tol : float
        Tolerance for convergence of the method.
    max_iter : int
        Maximum number of iterations.

    Methods
    -------
    anomaly_E(M: float, e: float, tol: float, max_iter: int) -> float
        Solves Kepler's equation: M = E - e*sin(E) to obtain the eccentric anomaly E.

    kepler(t: float, a: float, e: float, n: float) -> tuple
        Computes the position of the planet by solving Kepler's equation.

    """

    def __init__(self):
        self.M = 0
        self.e = 0
        self.E = 0
        self.a = 0
        self.r = None
        self.tol = 1e-8
        self.max_iter = 100
        self.n = 0
        self.time = None
        self.t0 = 0

    def anomaly_E(self, M:float) -> float:
            """
            Solves Kepler's equation: M = E - e*sin(E) to obtain the eccentric anomaly E
            for a given mean anomaly M and eccentricity e using the Newton-Raphson method.
            
            Parameters
            ----------
            M : float
                The mean anomaly in radians.
  
            Returns
            -------
            E : float
                The eccentric anomaly in radians.
            """
            
            from numpy import sin, cos, all, abs, copy
            # Initialize with E0 = M (good approximation for small e)
            E = copy(M)
            for i in range(self.max_iter):
                f = E - self.e * sin(E) - M  # ✅ Now correctly uses input M
                fp = 1 - self.e * cos(E)
                delta = -f / fp
                E += delta
                if all(abs(delta) < self.tol):
                    break
            return E
        
    def kepler(self) -> tuple:
            """
            Computes the position of the planet by solving Kepler's equation.
            
            Parameters
            ----------
  
            Returns
            -------
            r : array
                Cartesian coordinates of the planet in the orbital plane.
            """
            from numpy import arctan2, sqrt, sin, cos, zeros_like, array

            M = self.n * (self.time)
            E = self.anomaly_E(M)
            f = 2 * arctan2(sqrt(1 + self.e) * sin(E / 2), sqrt(1 - self.e) * cos(E / 2))
            r = self.a * (1 - self.e * cos(E))
            x = r * cos(f)
            y = r * sin(f)
            z = zeros_like(x)

            r = array([x, y, z]).T
            self.r = r

            return r
    
    def animate(self, save=False):
        """
        Smooth animation of the system motion with optimized performance.
        """
        from matplotlib import pyplot as plt
        from matplotlib.animation import FuncAnimation
        from IPython.display import HTML
        from matplotlib import rcParams

        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter3D(self.r[:, 0], self.r[:, 1], self.r[:, 2], 
                    c='black', marker='o', s=0.05, alpha=0.5)

        point, = ax.plot([], [], [], 'ro', markersize=10)
        ax.set_axis_off()
        ax.view_init(elev=90, azim=0)

        def update(frame):
            """Update function for animation"""
            x, y, z = self.r[frame, 0], self.r[frame, 1], self.r[frame, 2]
            point.set_data([x], [y])  # ✅ Ensure it's a sequence
            point.set_3d_properties([z])  # ✅ Ensure it's a sequence
            return point,

        frames = len(self.time)
        interval = 1000 / 30  # 🔥 30 FPS → 33ms per frame
        rcParams['animation.embed_limit'] = 2**128
        ani = FuncAnimation(fig, update, frames=frames, interval=interval, blit=True)

        if save:
            ani.save('media/kepler.gif', writer='imagemagick',fps=30)  # 🔥 Lower dpi for smaller size

        plt.close()
        return HTML(ani.to_jshtml())
 
       
    def load (self, path:str)->None:
         
        """
        Reads the parameters of the orbit from a txt file, and updates
        the attributes of the class

        Parameters
        ----------
        path : str
            The path to the txt file containing the parameters of the orbit

        Returns
        -------
        None
    
        """
        import pandas as pd
        with open(path, 'r') as file:
            lines = [line.strip() for line in file]
            data = dict()
            for line in lines:
                key, value = line.split('=')
                data[key] = float(value)
            self.a = data['a']
            self.b = data['b']
            self.e = data['e']
            self.c = data['c']
            self.p = data['p']
            self.q = data['q']
            self.n = data['n']
            self.f = data['f']
            self.E = data['E']
            self.M = data['M']
            self.T = data['T']
            self.t0 = data['t0']
            self.i = data['i']
            self.Ω = data['Ω']
            self.ω = data['ω']