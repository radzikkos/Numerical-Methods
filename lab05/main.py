from math import pow,sin,pi,fabs
import numpy as np

class MultiGridRelaxation():
    def __init__(self):
        self.delta = 0.2
        self.nx = 128
        self.ny = 128
        self.xMax = self.delta * self.nx
        self.yMax = self.delta * self.ny
        self.TOL = pow(10, -8)
        #self.V = np.zeros((self.nx + 1, self.ny + 1))
        self.k = [16, 8, 4, 2, 1]

    def stop_condition(self, V,k):
        S = 0.
        for i in range(0, self.nx - k + 1, k):
            for j in range(0, self.ny - k + 1, k):
                p1 = pow(k * self.delta, 2) / 2
                p2 = pow((V[i+k,j] - V[i,j])/(2 * k * self.delta) + (V[i+k, j+k] - V[i,j+k])/(2 * k * self.delta), 2)
                p3 = pow((V[i,j+k] - V[i,j])/(2 * k * self.delta) + (V[i+k, j+k] - V[i+k,j])/(2 * k * self.delta), 2)
                S += p1 * (p2 + p3)
        return S

    def edge_cases(self, V):
        for y in range(self.ny + 1):
            V[0,y] = sin(pi * self.delta * y / self.yMax)
        for x in range(self.nx + 1):
            V[x, self.ny] = (-1) * sin(2 * pi * self.delta * x / self.xMax)
        for y in range(self.ny + 1):
            V[self.nx, y] = sin(pi * self.delta * y / self.yMax)
        for x in range(self.nx + 1):
            V[x,0] = sin(2 * pi * self.delta * x / self.xMax)

    def condensing(self, V, k):
        for i in range(0, self.nx - k + 1, k):
            for j in range(0, self.ny - k + 1, k):
                V[i+k/2, j+k/2] = 1/4 * (V[i,j] + V[i+k,j] + V[i,j+k] + V[i+k, j+k])
                V[i+k, j+k/2] = 1/2 * (V[i+k, j] + V[i+k, j+k])
                V[i+k/2, j+k] = 1/2 * (V[i,j+k] + V[i+k,j+k])
                V[i+k/2, j] = 1/2 * (V[i,j] + V[i+k, j])
                V[i, j+k/2] = 1/2 * (V[i,j] + V[i, j+k])

    def discretization(self, V, k):
        for i in range(k, self.nx - k + 1, k):
            for j in range(k, self.ny - k + 1, k):
                V[i,j] = 1/4 * (V[i + k,j] + V[i - k,j] + V[i,j + k] + V[i,j - k])

    def solve(self):

        it = 0

        for step in self.k:
            sIt = np.zeros((1, 2))
            V = np.zeros((self.nx + 1, self.ny + 1))
            self.edge_cases(V)
            sIt[0,1] = self.stop_condition(V, step)

            while True:

                if(fabs( (sIt[0,1]- sIt[0,0]) / sIt[0,0] ) < self.TOL):
                    break
