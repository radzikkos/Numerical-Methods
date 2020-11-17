import math
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
class Relaxation:
    def __init__(self):
        self.epsilon = 1.
        self.delta = 0.1
        self.nx = 150
        self.ny = 100
        self.v1 = 10
        self.v2 = 0
        self.xMax = self.delta * self.nx
        self.yMax = self.delta * self.ny
        self.TOL = math.pow(10, -8)
        self.omegaG = [0.6, 1.0]
        self.sigmaX = 0.1 * self.xMax
        self.sigmaY = 0.1 * self.yMax

        #self.ro = np.zeros((self.nx, self.ny))
        self.omegaL = [1.0, 1.4, 1.8, 1.9]
    def ro1(self, x, y):
        part1 = - math.pow((x - 0.35 * self.xMax), 2) / math.pow(self.sigmaX, 2)
        part2 = - math.pow((y - 0.5 * self.yMax),2) / math.pow(self.sigmaY, 2)
        return math.exp(part1 + part2)

    def ro2(self, x, y):
        part1 = - math.pow((x - 0.65 * self.xMax), 2) / math.pow(self.sigmaX, 2)
        part2 = - math.pow((y - 0.5 * self.yMax), 2) / math.pow(self.sigmaY, 2)
        return -math.exp(part1 + part2)

    def fill_ro(self):
        tempRo = np.zeros((self.nx + 1, self.ny + 1))
        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                tempRo[i,j] = self.ro1(self.delta * i, self.delta * j) + self.ro2(self.delta * i, self.delta * j)
        return tempRo

    def stop_condition(self, V, ro):
        S = 0.
        for i in range(self.nx):
            for j in range(self.ny):
                S += math.pow(self.delta, 2) * (1/2 * (math.pow((V[i+1,j] - V[i,j]) / self.delta, 2)) + 1/2 * (math.pow((V[i,j+1]-V[i,j]) / self.delta, 2)) - ro[i,j] * V[i,j])
        return S

    def edge_cases(self, v):
        for j in range(self.nx + 1):
            v[j,0] = self.v1
            v[j,self.ny] = self.v2
        return v
    def global_relaxation_core(self, omega, vOld, vNew, ro):
        for i in range(1, self.nx):
            for j in range(1, self.ny):
                vNew[i,j] = 0.25 * ( vOld[i + 1,j] + vOld[i - 1,j] + vOld[i,j + 1] + vOld[i, j - 1] + math.pow(self.delta, 2) / self.epsilon * ro[i,j]);

        for j in range(1, self.ny):
            vNew[0,j] = vNew[1,j]
            vNew[self.nx, j] = vNew[self.nx - 1, j]

        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                vOld[i,j] = (1 - omega) * vOld[i,j] + omega * vNew[i,j]



    def err(self, vOld, i, j,ro):
        #ro = self.fill_ro()
        return ( (vOld[i+1,j] - 2*vOld[i,j] + vOld[i-1,j])/(math.pow(self.delta,2)) + (vOld[i,j+1] - 2*vOld[i,j] + vOld[i,j-1])/(math.pow(self.delta,2)) + ro[i,j]/self.epsilon );

    def global_relaxation(self):
        vOld = np.zeros((self.nx + 1, self.ny + 1))
        map = np.zeros((1, 3))  # to plot
        error = np.zeros((1, 3))  # to plot
        ro = self.fill_ro()
        for omega in self.omegaG:
            integral = np.zeros((1,2))  # to plot


            vNew = np.zeros((self.nx + 1, self.ny + 1))

            vNew = self.edge_cases(vNew)

            vOld = np.copy(vNew)
            sIt = np.zeros((1,2))
            sIt[0,1] = self.stop_condition(vNew, ro)

            it = 0
            while True:
                self.global_relaxation_core(omega, vOld, vNew, ro)
                sIt[0,0] = sIt[0,1]
                sIt[0,1] = self.stop_condition(vOld, ro)
                it += 1
                print(it)
                integral = np.vstack([integral, np.array([[it, sIt[0,0]]])])
                #print(integral, math.fabs( (sIt[0,0] - sIt[0,1])/ sIt[0,1]) - self.TOL)
                if( math.fabs( (sIt[0,0] - sIt[0,1])/ sIt[0,0]) < self.TOL ):
                    break
            plt.plot(integral[1:,0], integral[1:,1],label="omega="+(str)(omega))
            #print(vOld)
        plt.title("Integral for global relaxation")
        plt.xlabel("it")
        plt.ylabel("s(it)")
        plt.legend()
        plt.show()

        check = 0

        for i in range(self.nx + 1 ):
            for j in range(self.ny + 1):
                print(check)
                check += 1
                map = np.vstack([map, np.array([[self.delta * i, self.delta * j,vOld[i,j]]])])
                #print(i * j)
                if (i > 0 and j > 0 and i < self.nx and j < self.ny):
                    error = np.vstack([error, np.array([[self.delta * i, self.delta * j, self.err(vOld, i, j,ro)]])])
                else:
                    error = np.vstack([error, np.array([[self.delta * i, self.delta * j, 0.0]])])

        ax = plt.axes(projection='3d')

        ax.scatter3D(map[1:,0], map[1:,1], map[1:,2], cmap='viridis', linewidth=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z');
        plt.title("Map")
        plt.show()

        ax = plt.axes(projection='3d')
        ax.scatter3D(error[1:, 0], error[1:, 1], error[1:, 2], cmap='viridis', linewidth=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z');
        plt.title("Error")
        plt.show()

    def local_relaxation_core(self, omega, v, ro):
        for i in range(1,self.nx):
            for j in range(1,self.ny):
                v[i,j] = (1 - omega) * v[i,j] + (omega / 4.0) * (v[i + 1,j] + v[i - 1,j] + v[i,j + 1] + v[i,j - 1] + math.pow(self.delta, 2) / self.epsilon * ro[i,j]);

        for j in range(1,self.ny):
            v[0,j] = v[1,j]
            v[self.nx, j] = v[self.nx - 1, j]
        return v
    
    def local_relaxation(self):
        for omega in self.omegaL:
            integral = np.zeros((1,2))  # to plot
            sIt = np.zeros((1, 2))
            v = np.zeros((self.nx + 1, self.ny + 1))
            self.edge_cases(v)
            ro = self.fill_ro()
            sIt[0,1] = self.stop_condition(v, ro)

            it = 0

            while True:
                it += 1
                self.local_relaxation_core(omega, v, ro)
                sIt[0,0] = sIt[0,1]
                sIt[0,1] = self.stop_condition(v, ro)

                integral = np.vstack([integral, np.array([[it, sIt[0,0]]])])
                print(it)
                if (math.fabs((sIt[0, 0] - sIt[0, 1]) / sIt[0, 0]) < self.TOL):
                    break
            plt.plot(integral[1:, 0], integral[1:, 1], label="omega="+(str)(omega))
            # print(vOld)
        plt.title("Integral for local relaxation")
        plt.xlabel("it")
        plt.ylabel("s(it)")
        plt.legend()
        plt.show()
Relaxation().global_relaxation()
Relaxation().local_relaxation()


