'''
RRZ - kontrola kroku czasowego w problemach sztywnych.
'''
import math
import numpy as np
import matplotlib.pyplot as plt

x0 = 0.01
v0 = 0
deltaT0 = 1
tMax = 40
alpha = 5


class Methods():
    def __init__(self):
        self.TOL = (math.pow(10, -2), math.pow(10, -5))
        self.p = 2
        self.S = 0.75
    def rk2_method(self, xn, vn, deltaT, alpha):
        #Obylo sie bez jawnej definicji funckji f,g
        k1x = vn
        k1v = alpha * (1 - math.pow(xn, 2)) * vn - xn

        k2x = vn - deltaT * k1v
        k2v = alpha * (1 - math.pow(xn + deltaT * k1x, 2)) * (vn + deltaT * k1v) - (xn + deltaT * k1x)
        xnPlus1 = xn + deltaT / 2 * (k1x + k2x)
        vnPlus1 = vn + deltaT / 2 * (k1v + k2v)
        return xnPlus1, vnPlus1

    def trapezium_method(self, xn, vn, deltaT, alpha):
        def f(vn):
            return vn
        def g(xn, vn, alpha):
            return alpha * (1 - math.pow(xn,2)) * vn - xn

        def F(xnPlus1, xn, vnPlus1, vn, deltaT):
            return xnPlus1 - xn - deltaT / 2 * (f(vn) + f(vnPlus1))

        def G(xnPlus1, xn, vnPlus1, vn, deltaT, alpha):
            return vnPlus1 - vn - deltaT / 2 * (g(xn, vn, alpha) + g(xnPlus1, vnPlus1, alpha))

        def a21(deltaT, xnPlus1k, vnPlus1k, alpha):
            return (-1) * deltaT / 2 * ((-2) * alpha * xnPlus1k * vnPlus1k - 1)

        def a22(deltaT, xnPlus1k, alpha):
            return 1 - deltaT / 2 * alpha * (1 - math.pow(xnPlus1k, 2))
        sigma = math.pow(10, -10)
        xnPlus1 = xn
        vnPlus1 = vn
        a11 = 1
        a12 = -deltaT / 2.0
        while True:
            FResult = F(xnPlus1, xn, vnPlus1, vn, deltaT)
            GResult = G(xnPlus1, xn, vnPlus1, vn, deltaT, alpha)

            A21 = a21(deltaT, xnPlus1, vnPlus1, alpha)
            A22 = a22(deltaT, xnPlus1, alpha)

            denominator = (a11 * A22 - a12 * A21)

            deltaX = (-FResult * A22 - (-GResult) * a12) / denominator
            deltaV = (a11 * (-GResult) - A21 * (-FResult)) / denominator

            xnPlus1 += deltaX
            vnPlus1 += deltaV

            if math.fabs(deltaX) < sigma and math.fabs(deltaV) < sigma:
                break;
        return xnPlus1,vnPlus1

    #ta funkcja dziala poprawnie
    def time_step(self, deltaT0, x0, v0, tMax,alpha,fun,nameofMethod):
        dataTOL1 = data = np.array([[0, 0, 0, 0]])
        dataTOL2 =  np.array([[0, 0, 0, 0]])
        #Prowizoryczne zapisywanie do dwoch tablic, by pobrac z nich dane do wykresow
        i = 0
        for tol in self.TOL:
            i += 1
            t = 0
            deltaT = deltaT0
            xn = x0
            vn  = v0

            def ex(first, second):
                return (first - second) / (math.pow(2, self.p) - 1)
            def ev(first, second):
                return (first - second) / (math.pow(2, self.p) - 1)
            def max(first, second):
                if first >= second:
                    return  first
                return second
            j = 0
            while True:
            # stawiamy dwa kroki ∆t
                xvnPlus1 = fun(xn, vn, deltaT, alpha)
                xvnPlus2 = fun(xvnPlus1[0], xvnPlus1[1], deltaT, alpha)

            #stawiamy dwa kroki 2 * ∆t
                xvnPlus2v2 = fun(xn, vn, 2 * deltaT, alpha)

            #liczymy Ex, Ev
                Ex = ex(xvnPlus2[0], xvnPlus2v2[0])
                Ev = ev(xvnPlus2[1], xvnPlus2v2[1])

                maxValue = max(math.fabs(Ex), math.fabs(Ev))
                if maxValue < tol:
                    t = t + 2 * deltaT
                    xn = xvnPlus2[0]
                    vn = xvnPlus2[1]
                #zapus danych do pliku
                    if i == 1:
                        dataTOL1 = np.vstack([dataTOL1, np.array([[t, deltaT, xn, vn]])])
                    else:
                        dataTOL2= np.vstack([dataTOL2, np.array([[t, deltaT, xn, vn]])])
                #deltaT maleje, wiec czas przestaje sie zmieniac dla trapezow. Nie mam pojecia dlaczego tak sie dzieje
                deltaT = math.pow(self.S * tol / maxValue, 1 / (self.p + 1)) * deltaT
                if t >= tMax:
                    break
                #j += 1
                #print(deltaT)
                #if j > 10000:
                   # break
                #print(t , tMax)
        plt.plot(dataTOL1[1:, 0], dataTOL1[1:, 1], label="TOL =" + str(math.pow(10,-2)))
        plt.plot(dataTOL2[1:, 0], dataTOL2[1:, 1], label="TOL =" + str(math.pow(10,-5)))
        plt.title(nameofMethod)
        plt.xlabel("t")
        plt.ylabel("deltaT(t)")
        plt.legend()
        plt.show()

        plt.plot(dataTOL1[1:, 0], dataTOL1[1:, 2], label="TOL =" + str(math.pow(10, -2)))
        plt.plot(dataTOL2[1:, 0], dataTOL2[1:, 2], label="TOL =" + str(math.pow(10, -5)))
        plt.title(nameofMethod)
        plt.xlabel("t")
        plt.ylabel("x(t)")
        plt.legend()
        plt.show()

        plt.plot(dataTOL1[1:, 0], dataTOL1[1:, 3], label="TOL =" + str(math.pow(10, -2)))
        plt.plot(dataTOL2[1:, 0], dataTOL2[1:, 3], label="TOL =" + str(math.pow(10, -5)))
        plt.title(nameofMethod)
        plt.xlabel("t")
        plt.ylabel("v(t)")
        plt.legend()
        plt.show()

        plt.plot(dataTOL1[1:, 2], dataTOL1[1:, 3], label="TOL =" + str(math.pow(10, -2)))
        plt.plot(dataTOL2[1:, 2], dataTOL2[1:, 3], label="TOL =" + str(math.pow(10, -5)))
        plt.title(nameofMethod)
        plt.xlabel("x")
        plt.ylabel("v(x)")
        plt.legend()
        plt.show()

Obj = Methods()
#Rk2 dziala poprawnie
Obj.time_step(deltaT0, x0, v0, tMax, alpha,Obj.rk2_method,"RK2")
#Dziala
Obj.time_step(deltaT0, x0, v0, tMax, alpha,Obj.trapezium_method,"Trapezium")
