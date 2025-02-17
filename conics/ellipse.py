class Ellipse():
    """
    Class to represent an ellipse

    Attributes:
    a: float, major axis
    b: float, minor axis

    Methods:
    params: returns the parameters of the ellipse
    foci: returns the distance between the foci of the ellipse
    eccentricity: returns the eccentricity of the ellipse
    latus_rectum: returns the latus rectum of the ellipse
    directrix: returns the directrix of the ellipse
    area: returns the area of the ellipse
    perimeter: returns the perimeter of the ellipse
    Eq: returns the equation of the ellipse
    """

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def params(self)->tuple:
        """
        Parameters of the ellipse

        Returns:
        a: float
                major axis
        b: float
                minor axis

        c: float
                distance between the foci of the ellipse
        e: float
                eccentricity of the ellipse
        p: float
                semi-latus rectum of the ellipse
        q: float
                directrix of the ellipse
        """
        self.c = (self.a**2 - self.b**2)**0.5
        self.e = self.c / self.a
        self.p = self.a * (1 - self.e**2)
        self.q = self.b * (1 - self.e**2)
        return self.a, self.b, self.c, self.e, self.p, self.q

    def foci(self)->float:
        
        """Distance between the foci of the ellipse

        Returns:
        c: float
        """

        self.c = (self.a**2 - self.b**2)**0.5
        return self.c

    def eccentricity(self)->float:
        """
        Eccentricity of the ellipse
        
        Returns:
        e: float
        """
        self.c = self.foci()
        self.e = self.c / self.a
        return self.e

    def latus_rectum(self)->float:
        """
        Latus rectum of the ellipse

        Returns:
        p: float
        """
        self.c = self.foci()
        self.p = self.a * (1 - (self.c / self.a)**2)
        return self.p
    
    def directrix(self)->float:
        """
        Directrix of the ellipse

        Returns:
        q: float
        """
        self.c = self.foci()
        self.q = self.b * (1 - (self.c / self.a)**2)
        return self.q
    
    def area(self)->float:
        """
        Area of the ellipse

        Returns:
        area: float
        """
        self.area = 3.1416 * self.a * self.b
        return self.area
    
    def perimeter(self)->float:
        """
        Perimeter of the ellipse

        Returns:
        perimeter: float
        """
        self.perimeter = 2 * 3.1416 * ((self.a + self.b) / 2)
        return self.perimeter
    
    def plot(self)->None:
        """
        Plot the ellipse
        """
        import matplotlib.pyplot as plt
        import numpy as np
        t = np.linspace(0, 2 * np.pi, 100)
        x = self.a * np.cos(t)
        y = self.b * np.sin(t)
        plt.plot(x, y)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Ellipse')
        plt.grid()
        plt.show()

    def Eq(self):
        return f'x^2/{self.a**2} + y^2/{self.b**2} = 1'

    def __str__(self):
        return f'Ellipse with major axis {self.a} and minor axis {self.b}'
    
    def __repr__(self):
        return f'Ellipse({self.a}, {self.b})'
    






        


