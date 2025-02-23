class GCA():
    """
    The Guiding Centre Approximation (GCA) for approximanting the 2 body motion of a system

    Attributes:
    orb_params: Orbitals
                The parameters of the orbit
    semimajor (a) : float
                The semi-major axis of the orbit
    semiminor (b) : float
                The semi-minor axis of the orbit
    eccentricity (e) : float
                The eccentricity of the orbit
    focus (c) : float
                Distance between the foci of the orbit to the center
    latus_rectum (p) : float
                The semi-latus rectum of the orbit
    directrix (q) : float
                The directrix of the orbit
    mean_angular_velocity (n) : float
                The mean angular velocity of the orbit
    true_anomaly (f) : float
                The true anomaly of the orbit
    eccentric_anomaly (E) : float
                The eccentric anomaly of the orbit
    mean_anomaly (M) : float
                The mean anomaly of the orbit
    time_period (T) : float
                The time period of the orbit

    Methods:
    """
    from conics.orbitals import Orbitals
    from numpy import ndarray

    def __init__(self, orb_params=None)->None:
        from numpy import linspace, pi
        if orb_params is None:
            self.orb_params = self.Orbitals()
            self.a = None
            self.b = None
            self.e = None
            self.c = None
            self.p = None
            self.q = None
            self.n = None
            self.f = None
            self.E = None
            self.M = None
            self.time = linspace(0, 2*pi, 1000)
            self.M_sun = 1.989e30  # kg
            self.G = 6.67430e-11 # m^3/kg/s^2
            self.mu = self.G * self.M_sun # m^3/s^2
            self.AU = 1.496e11 # m
            self.rad2deg = 180/pi
            self.day2sec = 60*60*24 # s
            self.t0 = None
            self.inclination = None
            self.right_ascension = None
            self.argument_periapsis = None
            self.r_Gs = None
            self.v_Gs = None
            self.r_Ps = None
            self.v_Ps = None
            self.r = None
            self.v = None


        else:
            self.orb_params = orb_params
            self.a = orb_params.semimajor
            self.b = orb_params.semiminor
            self.e = orb_params.eccentricity
            self.c = orb_params.focus
            self.p = orb_params.latus_rectum
            self.q = orb_params.directrix
            self.n = orb_params.mean_angular_velocity
            self.f = orb_params.true_anomaly
            self.E = orb_params.eccentric_anomaly
            self.M = orb_params.mean_anomaly
            self.time = linspace(0, 2*pi, 1000)
            self.M_sun = 1.989e30  # kg
            self.G = 6.67430e-11 # m^3/kg/s^2
            self.mu = self.G * self.M_sun # m^3/s^2
            self.AU = 1.496e11 # m
            self.rad2deg = 180/pi
            self.day2sec = 60*60*24 # s
            self.t0 = orb_params.time_periapsis
            self.inclination = orb_params.inclination
            self.right_ascension = orb_params.right_ascension
            self.argument_periapsis = orb_params.argument_periapsis
            self.r_Gs = None
            self.v_Gs = None
            self.r_Ps = None
            self.v_Ps = None
            self.r = None
            self.v= None
            
    def proof(self)->None:
        """
        Proof of the GCA
        """
        from IPython.display import display, Markdown

        proof_text_1 = r"""Dado el caso de un problema de dos cuerpos gravitacional, si uno de estos orbita al otro siguiendo una trayectoria aproximadamente circular $(e << 1)$, se puede estudiar su movimiento descomponiendolo en una parte circular y otra eliptica descritas por dos puntos especiales."""
        
        proof_text_2 = r"""Un punto G en movimiento circular uniforme a una distancia (a) del centro, que corresponde con el semiejer mayor del cuerpo en su orbita, y una frecuencia angular media n resultado de la tercera ley de Kepler: $n = \sqrt{\frac{\mu}{a^3}}$"""
        
        proof_text_3 = r"""Ademas, se tiene un punto P orbitandolo con una trayectoria eliptica caracterizada por un semieje mayor y menor de longitud $2ae$ y $ae$ respectivamente, con (e) la excentricidad de la elipse. Este tambien orbita al punto G a la misma frecuencia angular media del mismo. El movimiento del cuerpo P, descrito desde un sistema de referencia centrado en G y a partir del diagrama en la figura 2.8 de MD, puede ser descrito en terminos de la anomalia verdadera y la anomlia media como: $x = r\cos{(f-M)} - a$ e $y = r\sin{(f-M)}$"""

        proof_text_4 = r"""Teniendo en cuenta la expasion del argumento de la funciones sinusoidales de cada coordenada, mostrado en la ecuacion 2.88 de MD, y el hecho de que por definicion de (e) cualquiera de sus potencias mayores a 1 tiende a cero, las expresiones para x e y se convierten a: $f - M \approx 2e\sin{M}$, $x \approx r\cos{(2e\sin{M})} - a$ e $y \approx r\sin{(2e\sin{M})}$"""

        proof_text_5 = r"""Estas expresiones pueden ser simplificadas realizando una nueva expansion en series de cosenos y senos para los terminos sinusoidales: $\cos{(2e\sin{M})} = 1 - 2e^{2}\sin{M}^2$ y $\sin{(2e\sin{M})} = 2e\sin{M}$"""
        
        proof_text_6 = r"""Ademas, de la ecuacion 2.20 se tiene que: $r = \frac{a(1-e^2)}{1 + e\cos{f}}$ el cual al expandir su denominador en torno a (e). Nuevamnete, siendo este un valor pequeño, se obtiene la expresion $(1-e\cos{f})$ en el denominador, mientras que en el numerador la potencia cuadratica desaparece por la misma razon, finalizando con la expresion: $r = \frac{a(1-e^2)}{1 + e\cos{f}} \approx a(1-e^2)\cdot(1-e\cos{f}) \approx a(1-e\cos{f})$"""

        proof_text_7 = r"""Y usando nuevamente el resultado de la ecuacion 2.88: $ f - M \approx 2e\sin{M}$ , se obtiene entonces que: $\cos{f} \approx \cos{M + 2e\sin{M}}$"""

        proof_text_8 = r"""Donde, dada la magnitud de e, el lado derecho adquiere la forma $\cos{M + \delta} \approx \cos{M} - \delta\sin{M}$, por lo cual: $\cos{f} \approx \cos{M} - \delta\sin{M} = \cos{M} - 2e\sin{M}\sin{M} = \cos{M} - 2e\sin{M}^2$"""

        proof_text_9 = r"""Teniendo el ultimo termino muy pequeño por el cuadrado del seno a la par de la multiplicacion por la excentricidad, la expresion final indica que, arpoximadamente, los cosenos de ambas anomalias es el mismo. $\cos{f} \approx \cos{M}$"""
        
        proof_text_10 = r"""Asi, reemplazando en las expresiones aproximadas de x e y, se obtiene:$x \approx a(1- e\cos{M})(1-2e^{2}\sin{M}^2) - a \approx - ae\cos{M}$ e $y \approx a(1- e\cos{M})2e\sin{M} \approx 2ae\sin{M}$"""

        proof_text_11 = r"""Finalmente, reemplazando en la ecuacion estandar de la elipse, se tiene que: $(\frac{- ae\cos{M}}{ae})^2  + (\frac{2ae\sin{M}}{2ae})^2 = \cos{M}^2 + \sin{M}^2 = 1$, cumpliendo con la identidad trigonometrica"""

        display(Markdown(proof_text_1))
        display(Markdown(proof_text_2))
        display(Markdown(proof_text_3))
        display(Markdown(proof_text_4))
        display(Markdown(proof_text_5))
        display(Markdown(proof_text_6))
        display(Markdown(proof_text_7))
        display(Markdown(proof_text_8))
        display(Markdown(proof_text_9))
        display(Markdown(proof_text_10))
        display(Markdown(proof_text_11))

    def decompose(self) -> tuple:
        """
        Decompose the motion of the system into two circular motions

        Parameters:
        time: np.ndarray
                The time array to decompose the motion

        Returns:
        x: np.ndarray
                The x component of the motion
        y: np.ndarray
                The y component of the motion
        """
        from numpy import sin, cos, array

        r_G = lambda a,n,t: array([a*cos(n*t), a*sin(n*t), 0])
        v_G = lambda a,n,t: array([-a*n*sin(n*t), a*n*cos(n*t), 0])

        r_P = lambda a,e, n, t: array([2*a*e*cos(n*t), a*e*sin(n*t), 0])
        v_P = lambda a,e, n, t: array([-a*n*e*sin(n*t), a*n*e*cos(n*t), 0])
        
        r_Gs = array([r_G(self.a,self.n,t) for t in self.time])
        v_Gs = array([v_G(self.a,self.n,t) for t in self.time])

        r_Ps = array([r_P(self.a,self.e, self.n,t) for t in self.time])
        v_Ps = array([v_P(self.a, self.e, self.n,t) for t in self.time])

        self.r_Gs = r_Gs
        self.v_Gs = v_Gs
        self.r_Ps = r_Ps
        self.v_Ps = v_Ps
        self.r = r_Gs + r_Ps
        self.v = v_Gs + v_Ps

        return r_Gs, v_Gs, r_Ps, v_Ps    

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

        ax.scatter3D(self.r_Gs[:, 0], self.r_Gs[:, 1], self.r_Gs[:, 2], 
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
            ani.save('media/comparison.gif', writer='imagemagick',fps=30)  # 🔥 Lower dpi for smaller size

        plt.close()
        return HTML(ani.to_jshtml())


    def __str__(self):
        return f"Orbitals: {self.orb_params}"
    
    def __params__(self):
        return dict(semimajor=self.a,
                    semiminor=self.b,
                    eccentricity=self.e,
                    focus=self.c,
                    latus_rectum=self.p,
                    directrix=self.q,
                    mean_angular_velocity=self.n,
                    true_anomaly=self.f,
                    eccentric_anomaly=self.E,
                    mean_anomaly=self.M
                    )
    