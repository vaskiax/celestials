import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pymcel as pc 
import plotly.graph_objects as go
import rebound as rb
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time
from astroquery.jplsbdb import SBDB

class System:
    def __init__(self, masses:dict, IC:dict, cannonics:dict):
        self.masses = masses
        self.IC = IC
        self.pos_pL = IC['p-L'][0]
        self.vel_pL = IC['p-L'][1]
        self.pos_epL = IC['e-pL'][0]
        self.vel_epL = IC['e-pL'][1]
        self.pos_sepL = IC['s-epL'][0]
        self.vel_sepL = IC['s-epL'][1]
        self.G = 6.67430e-11
        self.m_p = masses.loc['p'].values
        self.m_L = masses.loc['L'].values
        self.m_s = masses.loc['s'].values
        self.m_e = masses.loc['e'].values
        self.m_sepL = masses['mass'].sum()
        self.m_epL = self.m_sepL - self.m_s
        self.m_pL = self.m_epL - self.m_e
        self.M = np.array([self.m_s, self.m_e, self.m_p, self.m_L])
        self.UM = cannonics['UM']
        self.UL = cannonics['UL']
        self.UT = (self.UL**3/(self.G*self.UM))**0.5
        self.cannons = {'UM': self.UM, 'UL': self.UL, 'UT': self.UT}

        
        
    def solve(self, ts:np.ndarray) -> pd.DataFrame:
        """ 
        resuelve el problema de los n cuerpos en el tiempo ts usando la rutina ncuerpos_solucion
        de pymcel.

        Parameters
        ----------
        ts : np.ndarray
            Arreglo con los tiempos en los que se resolvera el problema de los n cuerpos.

        Returns
        -------
        pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos celestes en el tiempo.

        """
        from pymcel import ncuerpos_solucion

        rs, vs, rps, vps, quads = ncuerpos_solucion(self.IC, ts)

        df_rs = [pd.DataFrame(rs[i, :, :],
                                    columns = ['x', 'y', 'z']) for i in range(self.df.shape[0])]
        
        df_vs = [pd.DataFrame(vs[i, :, :],
                                    columns = ['vx', 'vy', 'vz']) for i in range(self.df.shape[0])]
        
        dfs = [pd.concat([df_rs[i], df_vs[i]], axis=1) for i in range(self.df.shape[0])]

        

        return pd.concat(dfs, axis=0, keys= [name for name in self.df.index])
    


def get_system_body_posvel(bodies:list[str], date:Time) -> pd.DataFrame:
    """
    El formato que recibe astropy es el siguiente:

    nombre_del_cuerpo = 'cuerpo_celeste'
    fecha = Time('AAAA-MM-DDTHH:MM:SS')

    La funcion get_body_barycentric_posvel recibe dos parametros, 
    el nombre de los cuerpos celestes y la fecha, y devuelve la posicion y
    velocidad del cuerpo celeste en el sistema baricentrico.

    Por lo tanto, esta funcion recibe una lista de cuerpos celestes y una fecha
    y devuelve un DataFrame con las posiciones y velocidades del sistema de cuerpos
    en la fecha requerida.

    Parameters
    ----------
    bodies : list[str]
        Lista con los nombres de los cuerpos celestes.
    date : Time
        Fecha en la que se quiere obtener la posicion y velocidad de los cuerpos celestes.

    Returns
    -------
    pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes.
        
    """

    pos_vel ={body: get_body_barycentric_posvel(body, date) for body in bodies}
    
    # Obtenemos las posiciones y velocidades de los cuerpos en forma vectorial
    pos_vel = {body: (pos_vel[body][0].xyz.to_value(),
                       pos_vel[body][1].xyz.to_value()) for body in pos_vel}
    
    return pd.DataFrame(pos_vel).transpose().rename(columns={0: 'position', 1: 'velocity'})



def get_solar_system_masses(bodies:list[str]) -> pd.DataFrame:
    """
    Obtiene las masas de los cuerpos mayores celestes del sistema solar en Kg

    Parameters
    ----------
    bodies : list[str]
        Lista con los nombres de los cuerpos celestes.

    Returns
    -------
    pd.DataFrame
        DataFrame con las masas de los cuerpos celestes en Kg.

    """

    masses = {
    'sun': 1.99e30,
    'mercury': 3.285e23,
    'venus': 4.867e24,
    'earth': 5.972e24,
    'moon': 7.342e22,
    'mars': 6.39e23,
    'jupiter': 1.898e27,
    'saturn': 5.683e26,
    'uranus': 8.681e25,
    'neptune': 1.024e26,
    }

    masses = pd.DataFrame({body: [masses[body]] for body in bodies})

    return masses.T.rename(columns={0: 'mass'})

    
def units_transform(df:pd.DataFrame) -> pd.DataFrame:
    """
    Transforma las unidades de las posiciones y velocidades de los cuerpos celestes
    de unidades astronomicas y dias a metros y segundos respectivamente.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes.
    
    Returns
    -------
    None

    """
    AU_TO_METERS = 149597870700  
    DAY_TO_SECONDS = 86400 

    df['position'] = df['position'].apply(lambda x: x * AU_TO_METERS)
    df['velocity'] = df['velocity'].apply(lambda x: x * AU_TO_METERS / DAY_TO_SECONDS)

    return df



def evolution(df_system:pd.DataFrame, masses:dict, times:np.ndarray, G:float = 6.67430e-11) -> np.ndarray:

    """
    Calcula la evolucion de las posiciones de los cuerpos celestes en el tiempo.

    Parameters
    ----------
    df_system : pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes.
    masses : dict
        Diccionario con las masas de los cuerpos celestes.
    times : np.ndarray
        Arreglo con los tiempos en los que se calculara la evolucion de las posiciones.
    
    Returns
    -------
    np.ndarray
        Arreglo con las posiciones de los cuerpos celestes en el tiempo.

    """
    sim = rb.Simulation()
    sim.G = G
    #sim.dt = 86400
    for i in df_system.index:
        sim.add(m = masses.loc[f'{i}', 'mass'],
                x = df_system.loc[i, 'x'],
                y = df_system.loc[i, 'y'],
                z = df_system.loc[i, 'z'],
                vx = df_system.loc[i, 'vx'],
                vy = df_system.loc[i, 'vy'],
                vz = df_system.loc[i, 'vz'])
        
    num_particles = len(sim.particles)
    positions = np.zeros((len(times), num_particles, 3))
    velocities = np.zeros((len(times), num_particles, 3))

    for i, time in enumerate(times):
        sim.integrate(time)
        for j in range(num_particles):
            positions[i, j] = [sim.particles[j].x, sim.particles[j].y, sim.particles[j].z]
            velocities[i, j] = [sim.particles[j].vx, sim.particles[j].vy, sim.particles[j].vz]

    return positions, velocities



def plot_evolution(positions:np.ndarray, body_names:list, typ = 'a') -> None:
    """ 
    
    Grafica la evolucion de las posiciones de los cuerpos celestes en el tiempo.

    Parameters
    ----------
    positions : np.ndarray
        Arreglo con las posiciones de los cuerpos celestes en el tiempo.
    body_names : list
        Lista con los nombres de los cuerpos celestes.

    Returns
    -------
    None

    """
    if typ == 'a':

        fig = go.Figure()
        num_particles = positions.shape[1]
        for i in range(num_particles):
            fig.add_trace(go.Scatter3d(x = positions[:, i, 0],
                                        y = positions[:, i, 1],
                                        z = positions[:, i, 2],
                                            mode = 'lines',
                                            name = body_names[i]))
        fig.show()

    if typ == 'z':
             
            fig = go.Figure()
            num_particles = positions.shape[0]
            for i in range(num_particles):
                fig.add_trace(go.Scatter3d(x = positions[i, :, 0],
                                            y = positions[i, :, 1],
                                            z = positions[i, :, 2],
                                                mode = 'lines',
                                                name = body_names[i]))
                

            fig.update_layout(
                scene=dict(
                    xaxis=dict(
                        backgroundcolor="black",
                        gridcolor="black",
                        showbackground=True,
                        zerolinecolor="black",
                        color="white",
                        showticklabels=False,
                        title_text = ''
                    ),
                    yaxis=dict(
                        backgroundcolor="black",
                        gridcolor="black",
                        showbackground=True,
                        zerolinecolor="black",
                        color="white",
                        showticklabels=False,
                        title_text = ''
                    ),
                    zaxis=dict(
                        backgroundcolor="black",
                        gridcolor="black",
                        showbackground=True,
                        zerolinecolor="black",
                        color="white", 
                        showticklabels=False,
                        title_text = ''
                    )
                ),
                paper_bgcolor="black",
                plot_bgcolor="black"
            )

            fig.show()


def plot_evolution2D (positions:np.ndarray, body_names:list, typ = 'a') -> None:
    """ 
    
    Grafica la evolucion de las posiciones de los cuerpos celestes en el tiempo en 2D.

    Parameters
    ----------
    positions : np.ndarray
        Arreglo con las posiciones de los cuerpos celestes en el tiempo.
    body_names : list
        Lista con los nombres de los cuerpos celestes.

    Returns
    -------
    None

    """
    if typ == 'a':
        fig = plt.figure(figsize=(10,8)) 
        num_particles = positions.shape[1]
        for i in range(num_particles):
            plt.plot(positions[:, i, 0], positions[:, i, 1], label = body_names[i])
        plt.legend()
        plt.axis('equal')
        plt.show()

    if typ == 'z':
        fig = plt.figure(figsize=(10,8)) 
        num_particles = positions.shape[0]
        for i in range(num_particles):
            plt.plot(positions[i, :, 0], positions[i, :, 1], label = body_names[i])
        plt.legend()
        plt.axis('equal')
        plt.show()



def separate_df(df:pd.DataFrame) -> pd.DataFrame:
    """
    Separa las posiciones y velocidades de los cuerpos celestes en formato vectorial
    a formato por componente.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes en formato vectorial.
    
    Returns
    -------
    pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes por componente.

    """

    rs = np.stack(df['position'].values)
    vs = np.stack(df['velocity'].values)

    dfs = pd.concat([pd.DataFrame({
                                      'x': rs[:, i],
                                      'y': rs[:, i+1],
                                      'z': rs[:, i+2],
                                      'vx': vs[:, i],
                                      'vy': vs[:, i+1],
                                      'vz': vs[:, i+2]})

                                      for i in range(0, len(rs[0]), 3)],
                                      axis=1).set_index(df.index)
    return dfs


def aggregate (df:pd.DataFrame) -> pd.DataFrame:
    """
    Agrega las posiciones y velocidades de los cuerpos celestes en formato por componente
    a formato vectorial.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes por componente.
    
    Returns
    -------
    pd.DataFrame
        DataFrame con las posiciones y velocidades de los cuerpos celestes en formato vectorial.

    """

    rs = np.stack([df['x'].values, df['y'].values, df['z'].values], axis=1)
    vs = np.stack([df['vx'].values, df['vy'].values, df['vz'].values], axis=1)

    dfs = pd.DataFrame({'position': rs.tolist(),
                        'velocity': vs.tolist()}, index = df.index)
    return dfs



def L(df:pd.DataFrame, masses:pd.DataFrame) -> np.ndarray:

        """
        Calcula el momento angular de cada cuerpo en cada instante de tiempo.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos.
        masses : pd.DataFrame
                DataFrame con las masas de los cuerpos.

        Returns
        -------

        np.ndarray
                Momento angular del sistema en cada instante de tiempo.

        """
        
        if len(df.index.get_level_values(0).unique()) == df[df.columns[0]].size:
            
            l = np.array([np.cross(df.loc[body, 'position'],
                             df.loc[body, 'velocity'])*masses.loc[body, 'mass']
                             for body in masses.index
                             
                             ]).sum(axis=0)
        else:

            l = [np.cross(df.loc[body, 'position'].iloc[i],
                                df.loc[body, 'velocity'].iloc[i])*masses.loc[body, 'mass']
                                
                                for body in masses.index
                                for i in range(len(df.loc[body]))]


        return l

def RCM(df:pd.DataFrame, masses:pd.DataFrame, PCM_:np.ndarray) -> np.ndarray:

        """
        Calcula la cuadratura del centro de masa de cada cuerpo en cada instante de tiempo.
        
        Sobre un sistema con un un solo registro de tiempo, devuelve el valor total del sistema,
        mientras que para un sistema de varios registros, devuelve 
        el centro de masa de cada particula del sistema en cada tiempo 

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos.
        masses : pd.DataFrame
                DataFrame con las masas de los cuerpos.
        PCM_=0 : np.ndarray
                vector momento lineal del sistema

        Returns
        -------

        np.ndarray
                Momento lineal del sistema en cada instante de tiempo.

        """
        

        if len(df.index.get_level_values(0).unique()) == df[df.columns[0]].size:
            
            p = [np.array(np.array(df.loc[body, 'position'])*masses.loc[body,'mass'] - 
                                    PCM_.loc[body]) / masses.loc['mass'].sum()
                                    
                    for body in masses.index]
        else:
            p = [np.array(np.array(df.loc[body, 'position'].iloc[time])*masses.loc[body,'mass'] - 
                                    PCM_.loc[body, time]) / masses['mass'].sum()

                    for body in masses.index
                    for time in range(len(df.loc[body]))]

        return p

def K(df:pd.DataFrame, masses:pd.DataFrame) -> np.ndarray:

        """
        Calcula la energia cinetica de cada cuerpo en cada instante de tiempo.
        
        Sobre un sistema con un un solo registro de tiempo, devuelve el momento
        lineal total, mientras que para un sistema de varios registros, devuelve 
        la energia cinetica individual de cada particula del sistema en cada tiempo 

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos.
        masses : pd.DataFrame
                DataFrame con las masas de los cuerpos.

        Returns
        -------

        np.ndarray
                Energia cinetica del sistema en cada instante de tiempo.

        """
        

        if len(df.index.get_level_values(0).unique()) == df[df.columns[0]].size:
            p = np.array([0.5*masses.loc[body,'mass']*np.linalg.norm(df.loc[body, 'velocity'])**2   
                        for body in masses.index]).sum(axis=0)
        else:
            p = [0.5*masses.loc[body,'mass']*np.linalg.norm(df.loc[body, 'velocity'].iloc[time])**2   
                        for body in masses.index for time in range(len(df.loc[body]))]

        return p


def U(df:pd.DataFrame, masses:pd.DataFrame, G = 1) -> np.ndarray:

        """
        Calcula la energia potencial de cada cuerpo en cada instante de tiempo.
        
        Sobre un sistema con un un solo registro de tiempo, devuelve el momento
        lineal total, mientras que para un sistema de varios registros, devuelve 
        la energia potencial individual de cada particula del sistema en cada tiempo 

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos.
        masses : pd.DataFrame
                DataFrame con las masas de los cuerpos.

        Returns
        -------

        np.ndarray
                Energia potencial del sistema en cada instante de tiempo.

        """
        
        df_bodies = lambda df: df.index.get_level_values(0).unique()

        if len(df.index.get_level_values(0).unique()) == df[df.columns[0]].size:
            mm = masses['mass'].values
            rr = np.stack([df.loc[body, 'position'] for body in masses.index], axis=1)

            p = np.array([sum([mm[i]*mm[j]/np.linalg.norm(rr[k][i]-rr[k][j])
                                
                                for k in range(len(rr))
                                for j in range(len(mm)) if i != j])
                                for i in range(len(mm))]).sum()
        else:
            
            time_r = len(df.loc[df_bodies(df)[0]])
            p = np.stack([np.array(
                            [
                                -G *
                                masses.loc[body_a, 'mass'] * masses.loc[body_b, 'mass'] /
                                np.linalg.norm(
                                    np.array(df.loc[(body_a, time), 'position']) -
                                    np.array(df.loc[(body_b, time), 'position'])
                                ) 
                                
                                for body_a in masses.index
                                for body_b in masses.index.drop(body_a)
                                
                            
                            ]
                
                        ).reshape(len(df_bodies(df)), 3).sum(axis=1) / 2
                        
                        for time in range(time_r)])
            
            p = np.concatenate([p[:, i] for i in range(len(p[0]))])


        return p


def V(Us:np.ndarray, Ks:np.ndarray) -> np.ndarray:

        """
        Calcula el virial del sistema en cada instante de tiempo.

        Parameters
        ----------
        Us : np.ndarray
            Array con las energias potenciales del sistema en cada tiempo.
        Ks : np.ndarray
            Array con las energias cineticas del sistema en cada tiempo.

        Returns
        -------

        np.ndarray
                Virial del sistema en cada instante de tiempo.

        """
        

        if type(Us) == float or type(Us) == np.float64:
            p = Us + 2*Ks
        else:
            p = [u + 2*k for u,k in zip(Us, Ks)]

        return p


def P(df:pd.DataFrame, masses:pd.DataFrame) -> np.ndarray:

        """
        Calcula el momento lineal de cada cuerpo en cada instante de tiempo.
        
        Sobre un sistema con un un solo registro de tiempo, devuelve el momento
        lineal total, mientras que para un sistema de varios registros, devuelve 
        el momento lineal individual de cada particula del sistema en cada tiempo 

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos.
        masses : pd.DataFrame
                DataFrame con las masas de los cuerpos.

        Returns
        -------

        np.ndarray
                Momento lineal del sistema en cada instante de tiempo.

        """
        

        if len(df.index.get_level_values(0).unique()) == df[df.columns[0]].size:
            p = np.array([np.array([

                        np.array(df.loc[body, 'velocity'])*masses.loc[body, 'mass']]).sum(axis=0)
                        
                        for body in df.index.get_level_values(0).unique()]).sum(axis=0)
        
        else:
            p = [np.array([

                        np.array(df.loc[(body, time), 'velocity'])*masses.loc[body, 'mass']]).sum(axis=0)
                        
                        for body in df.index.get_level_values(0).unique()
                        for time in range(len(df.loc[df.index.get_level_values(0).unique()[0]]))

                        ]

        return p


def cube2df (cube_rs:np.ndarray, cube_vs:np.ndarray, bodies:list, typ = 'a') -> pd.DataFrame:

    """ 
    Función que convierte un cubo de posiciones o velocidades en un DataFrame

    Parameters
    ----------
    cube : np.ndarray
        Cubo de posiciones o velocidades
    bodies : list
        Lista de cuerpos

    Returns
    -------
    pd.DataFrame
        DataFrame con las posiciones o velocidades de los cuerpos
        
    """
    if typ == 'a':
        N =  cube_rs.shape[1]
        dfs_rs = pd.concat([pd.DataFrame(cube_rs[:, i, :],
                            columns = ['x', 'y', 'z']) for i in range(N)],
                            axis=0, keys=[name for name in bodies])
        
        dfs_vs = pd.concat([pd.DataFrame(cube_vs[:, i, :],
                            columns = ['vx', 'vy', 'vz']) for i in range(N)],
                            axis=0, keys=[name for name in bodies])
        
        return pd.concat([dfs_rs, dfs_vs], axis=1)
        
    if typ == 'z':
        N =  cube_rs.shape[0]
        dfs_rs = pd.concat([pd.DataFrame(cube_rs[i, :, :],
                            columns = ['x', 'y', 'z']) for i in range(N)],
                            axis=0, keys=[name for name in bodies])
        
        dfs_vs = pd.concat([pd.DataFrame(cube_vs[i, :, :],
                            columns = ['vx', 'vy', 'vz']) for i in range(N)],
                            axis=0, keys=[name for name in bodies])

        return pd.concat([dfs_rs, dfs_vs], axis=1)
        
    if typ == 'c':
        N =  cube_rs.shape[0]
        dfs_rs = pd.concat([pd.DataFrame(cube_rs[i, :, :],
                            columns = ['x', 'y', 'z']) for i in range(N)],
                            axis=0, keys=[name for name in bodies])


        return dfs_rs
    

def Quad(df:pd.DataFrame, param:str) -> np.ndarray:
    
    """
        Calcula la cuadratura de un sistema de cuerpos en todos los tiempos
        según el parámetro dado.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame con las posiciones y velocidades de los cuerpos.
        param : str
            Parámetro a evaluar.

        Returns
        -------
        np.ndarray
            Array con las cuadraturas del sistema para cada tiempo

    """


    if param in ['PCM', 'RCM', 'L']:

        temp = np.array([np.array(df.loc[body, param].iloc[time])
                        for time in range(len(df.loc[df.index.get_level_values(0).unique()[0]]))
                        for body in df.index.get_level_values(0).unique()])

        temp = temp.reshape(len(df.loc[df.index.get_level_values(0).unique()[0]]),
                            len(df.index.get_level_values(0).unique()), 3)

    elif param in ['K', 'U']:

        temp = np.array([np.array(df.loc[body, param].iloc[time])
                        for time in range(len(df.loc[df.index.get_level_values(0).unique()[0]]))
                        for body in df.index.get_level_values(0).unique()])
        
        temp = temp.reshape(len(df.loc[df.index.get_level_values(0).unique()[0]]),
                     len(df.index.get_level_values(0).unique()))
    
    return [t.sum(axis=0) for t in temp]

years = lambda x: x*365*24*3600

sistema_z = lambda df, masses, bodies: [dict(m=masses.loc[f'{i}'].values,
                                             r = df.loc[i, ['x', 'y', 'z']].values,
                                                v = df.loc[i, ['vx', 'vy', 'vz']].values) for i in bodies]


def Orbitals (astro:str, params:list[str]) -> pd.DataFrame:

    """
    Función que obtiene los elementos orbitales de un cuerpo celeste.

    Parameters
    ----------
    astro : str
        Nombre del cuerpo celeste.
    params : list
        Lista con los parametros a obtener.

    Returns
    -------
    pd.DataFrame
        DataFrame con los elementos orbitales del cuerpo celeste.

    """

    sbdb = pd.DataFrame(SBDB.query(astro))

    return pd.Series({f'{name}': sbdb.loc['elements', 'orbit'][name].value
                      if name != 'e' else sbdb.loc['elements', 'orbit'][name]
                       for name in params})


def rot_matrix(theta:float, simetry:str) -> np.ndarray:

    """
    Función que obtiene la matriz de rotación en 3D.

    Parameters
    ----------
    theta : float
        Angulo de rotación.
    simetry : str
        Eje de simetría.

    Returns
    -------
    np.ndarray
        Matriz de rotación en 3D.

    """


    z_ = (0,1,3,4)
    y_ = (0,2,6,8)
    x_ = (4,5,7,8)

    kernel = np.eye(3).reshape(9)
    nucleoids = (np.cos(theta), np.sin(theta), -np.sin(theta), np.cos(theta))

    if simetry == 'z':
        kernel[[z_]] = nucleoids
        return kernel.reshape(3,3)
    
    if simetry == 'y':
        kernel[[y_]] = nucleoids
        return kernel.reshape(3,3)
    
    if simetry == 'x':
        kernel[[x_]] = nucleoids
        return kernel.reshape(3,3)
