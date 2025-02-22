class Orbitals():
    def __init__(self):
        self.semimajor = 0
        self.semiminor = 0
        self.eccentricity = 0
        self.focus = 0
        self.latus_rectum = 0
        self.directrix = 0
        self.mean_angular_velocity = 0
        self.true_anomaly = 0
        self.eccentric_anomaly = 0
        self.mean_anomaly = 0
        self.time_period = 0
        self.time_periapsis = 0
        self.inclination = 0
        self.right_ascension = 0
        self.argument_periapsis = 0


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
            self.semimajor = data['a']
            self.semiminor = data['b']
            self.eccentricity = data['e']
            self.focus = data['c']
            self.latus_rectum = data['p']
            self.directrix = data['q']
            self.mean_angular_velocity = data['n']
            self.true_anomaly = data['f']
            self.eccentric_anomaly = data['E']
            self.mean_anomaly = data['M']
            self.time_period = data['T']
            self.time_periapsis = data['t0']
            self.inclination = data['i']
            self.right_ascension = data['Ω']
            self.argument_periapsis = data['ω']
            



        
        