import numpy as np
from numpy import sin, cos
class Tensor_form_scara():
    def __init__(self, ac1, a1, ac2, a2, m1, i1zz, m2, i2zz, m3, i3zz, th1, th2, thv1, thv2):
        self.ac1 = ac1
        self.a1 = a1
        self.ac2 = ac2
        self.a2 = a2
        self.m1 = m1
        self.i1zz = i1zz
        self.m2 = m2
        self.i2zz = i2zz
        self.m3 = m3
        self.i3zz = i3zz
        self.th1 = np.deg2rad(th1)
        self.th2 = np.deg2rad(th2)
        self.thv1 = np.deg2rad(thv1)
        self.thv2 = np.deg2rad(thv2)

    def mass_matrix(self):
        m11 = (self.m1*self.ac1**2)+(self.m2+self.m3)*self.a2+(2*self.m3*self.a1*self.a2*cos(self.th2))+(self.m2*self.ac2**2)+(self.m3*self.a2**2)+(
                2*self.m2*self.a1*self.ac2*cos(self.th2))+(2*self.a1*cos(self.th2*(self.m2*self.ac2+self.m3*self.a2)))+self.i1zz+self.i2zz+self.i3zz
        m12 = (self.m2*self.ac2**2)+(self.m2*self.a1*self.ac2*cos(self.th2))+(self.m3*self.a2**2)+(self.m3*self.a1*self.a2*cos(self.th2))+(
                cos(self.th2)*self.a1*(self.m2*self.ac2+self.m3*self.a2))+self.i2zz+self.i3zz
        m13 = 0
        m14 = self.i3zz
        m21 = self.m2*self.ac2**2+(self.m2*cos(self.th2)*self.a1*self.ac2)+(self.m3*self.a2**2)+(self.m3*cos(self.th2)*self.a1*self.a2)+self.i3zz
        m22 = (self.m2*self.ac2**2)+(self.m3*self.a2**2)+self.i3zz
        m23 = 0
        m24 = self.i3zz
        m31 = 0
        m32 = 0
        m33 = self.m3
        m34 = 0
        m41 = self.i3zz
        m42 = self.i3zz
        m43 = 0
        m44 = self.i3zz

        M = np.array([[m11, m12, m13, m14],
            [m21, m22, m23, m24],
            [m31, m32, m33, m34],
            [m41, m42, m43, m44]])
        #print(f'{m11}   {m12}   {m13}   {m14}\n{m21}    {m22}   {m23}   {m24}\n{m31}    {m32}   {m33}   {m34}\n{m41}    {m42}   {m43}   {m44}')
        print(f'MASS MATRIX:\n {M}')
    
    def kinetics_matrix(self):
        c11 = (-2*self.m2*sin(self.th1)*self.a1*self.ac2)-(2*self.m3*sin(self.th2)*self.a1*self.a2)
        c12 = -(self.m2*sin(self.th2)*self.a1*self.ac2)-(self.m3*sin(self.th2)*self.a1*self.a2)
        c13 = 0
        c14 = 0
        c21 = -(self.m2*self.thv2*self.a1*self.ac2)-(sin(self.th2)*self.a1*((self.m3*self.thv2*self.ac2)+(self.m3*self.thv1*self.a2)))
        c22 = (self.a1*sin(self.th2))*((self.m2*self.ac2*self.thv1)+(self.m3*self.a2*self.thv1))
        
        C = np.array([[c11, c12, c13, c14],
            [c21, c22, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
        print(f'KINETICS MATRIX:\n {C}')
    
    def gravities_matrix(self):
        G = np.array([[0],
            [0],
            [-self.m3*9.81],
            [0]])
        print(f'GRAVITIES MATRIX:\n {G}')

tensor = Tensor_form_scara(0.15, 0.3, 0.1071, 0.251, 1.865861, 0.026426, 6.438328, 0.065753, 0.147655, 0.000011, 30, 30, 10, 15)
print(tensor.mass_matrix())
print(tensor.kinetics_matrix())
print(tensor.gravities_matrix())