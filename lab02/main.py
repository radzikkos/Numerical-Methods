import math
import numpy as np
import matplotlib.pyplot as plt

class TrapeziumMethod():
    def __init__(self):
        self.betha = 0.001
        self.N = 500
        self.gamma = 0.1
        self.tMax = 100
        self.deltaT = 0.1
        self.u0 = 1
        self.TOL = math.pow(10,-6)
        self.maxIteration = 20

    def f(self, un):
        alfa = self.betha * self.N - self.gamma
        return alfa * un - self.betha * math.pow(un,2)

    def df(self,un):
        alfa = self.betha * self.N - self.gamma
        return alfa - 2 * self.betha * un

    def pickard_method(self):

        un = self.u0
        unPlusOne = self.u0
        unArray = np.array([[0,self.u0]]) #time and un
        actualT = self.deltaT

        while actualT <= self.tMax:
            actualMi = 0
            un = unPlusOne
            uMi = 0.
            while math.fabs(unPlusOne - uMi) > self.TOL and actualMi < self.maxIteration :
                uMi = unPlusOne
                unPlusOne = un + (self.deltaT / 2.) * (self.f(un) * self.f(uMi) )
                actualMi += 1
            unArray = np.vstack([unArray, np.array([[actualT, unPlusOne]])])

            actualT += self.deltaT
        plt.plot(unArray[:, 0], unArray[:, 1], label="u(t)")
        plt.plot(unArray[:, 0], self.N - unArray[:, 1],label="z(t)")
        plt.title("Pickard Method")
        plt.xlabel("t")
        plt.ylabel("u(t) , z(t)")
        plt.legend()
        plt.show()

    def newton_method(self):
        un = self.u0
        unPlusOne = self.u0
        unArray = np.array([[0, self.u0]])  # time and un
        actualT = self.deltaT

        while actualT <= self.tMax:
            actualMi = 0
            un = unPlusOne
            uMi = 0.
            while math.fabs(unPlusOne - uMi) > self.TOL and actualMi < self.maxIteration:
                uMi = unPlusOne
                unPlusOne = uMi - (uMi - un - self.deltaT / 2 * (self.f(un) * self.f(uMi)) / (1 - self.deltaT / 2 * self.df(uMi)))
                actualMi += 1
            unArray = np.vstack([unArray, np.array([[actualT, unPlusOne]])])

            actualT += self.deltaT
        plt.plot(unArray[:, 0], unArray[:, 1], label="u(t)")
        plt.plot(unArray[:, 0], self.N - unArray[:, 1], label="z(t)")
        plt.title("Newton Method")
        plt.xlabel("t")
        plt.ylabel("u(t) , z(t)")
        plt.legend()
        plt.show()


class RK2:
    def __init__(self):
        self.betha = 0.001
        self.N = 500
        self.gamma = 0.1
        self.tMax = 100
        self.deltaT = 0.1
        self.u0 = 1
        self.TOL = math.pow(10,-6)
        self.maxIteration = 20
        self.alpha = self.betha * self.N - self.gamma

    def f1(self, U1, U2, un, a11, a12):
        return U1 - un - self.deltaT * (a11 * (self.alpha * U1 - self.betha * math.pow(U1,2)) + a12 * (self.alpha * U2 - self.betha * math.pow(U2,2)))

    def f2(self, U1, U2, un, a21, a22):
        return U2 - un - self.deltaT * (a21 * (self.alpha * U1 - self.betha * math.pow(U1,2)) + a22 * (self.alpha * U2 - self.betha * math.pow(U2,2)))

    def m11(self, a11, U1):
        return 1 - self.deltaT * a11 * (self.alpha - 2 * self.betha * U1)

    def m12(self, a12, U2):
        return (-1) * self.deltaT * a12 * (self.alpha - 2 * self.betha * U2)

    def m21(self, a21, U1):
        return (-1) * self.deltaT * a21 * (self.alpha - 2 * self.betha * U1)

    def m22(self, a22, U2):
        return 1 - self.deltaT * a22 * (self.alpha - 2 * self.betha * U2)

    def deltaU1(self, F1, F2, m22, m11, m12, m21):
        return ( F2 * m12 - F1 * m22 ) / ( m11 * m22 - m12 * m21)

    def deltaU2(self, F1, F2, m22, m11, m12, m21):
        return ( F1 * m21 - F2 * m11 ) / (m11 * m22 - m12 * m21)

    def f(self, un):
        return self.alpha * un - self.betha * math.pow(un, 2)

    def solve(self):
        a22 = a11 = 1/4
        a12 = 1/4 - math.sqrt(3) / 6
        a21 = 1/4 + math.sqrt(3) / 6
        b = 1/2
        unArray = np.array([[0, self.u0]])  # time and un
        actualT = self.deltaT
        unPlusOne = self.u0

        while actualT <= self.tMax:
            un = unPlusOne
            U1Mi = 0
            U2Mi = 0
            actualMi = 0
            U1 = un
            U2 = un
            while math.fabs(U1 - U1Mi) > self.TOL or math.fabs(U2 - U2Mi) > self.TOL and actualMi < self.tMax:
                U1Mi = U1
                U2Mi = U2
                deltaTU1 = self.deltaU1(self.f1(U1, U2, un, a11, a12), self.f2(U1, U2, un, a21, a22),
                                        self.m22(a22, U2), self.m11(11,U1), self.m12(a12,U2), self.m21(a21,U1))
                deltaTU2 = self.deltaU2(self.f1(U1, U2, un, a11, a12), self.f2(U1, U2, un, a21, a22),
                                        self.m22(a22, U2), self.m11(11,U1), self.m12(a12,U2), self.m21(a21,U1))
                U1 = U1Mi + deltaTU1
                U2 = U2Mi + deltaTU2
                actualMi += 1
            unPlusOne = un + self.deltaT * (b * self.f(U1)  + b * self.f(U2))
            unArray = np.vstack([unArray, np.array([[actualT, unPlusOne]])])
            actualT += self.deltaT
        plt.plot(unArray[:,0], unArray[:,1], label="u(t)")
        plt.plot(unArray[:, 0],self.N - unArray[:, 1], label="z(t)")
        plt.title("RK2 Method")
        plt.xlabel("t")
        plt.ylabel("u(t) , z(t)")
        plt.legend()
        plt.show()
TrapeziumMethod().pickard_method()
TrapeziumMethod().newton_method()
RK2().solve()