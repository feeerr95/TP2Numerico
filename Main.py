import math
from matplotlib import pyplot as plt

"""Posicion de la Tierra en un sistema X-Y fijo con unidades de metros"""
X1 = -4.670 * 10 ** 6
Y1 = 0
"""Posicion de la Luna en un sistema X-Y fijo con unidades de metros"""
X2 = 379.7 * 10 ** 6
Y2 = 0
'Radio de la Tierra y de la Luna respectivamente'
R1 = 6.731 * 10 ** 6
R2 = 1.738 * 10 ** 6
'Masa de la Tierra y de la Luna respectivamente'
M1 = 5972 * 10 ** 21
M2 = 73.48*10**21
'Constante universal de gravitacion'
G = 6.674 * 10 ** (-11)
'Distancia inicial entre la sonda y la superficie de la Tierra'
H0 = 0.373 * 10 ** 6
'Velocidad angular de giro del sistema Tierra-Luna alrededor de su CM'
W = 4.236*10**(-7)

Vx = []
Vy = []
X = []
Y = []
H = 1
#N = int(((2 * math.pi * (R1+H0)) / math.sqrt((G*M1)/(R1+H0))) / H)
N = 1200000

"""Calcula la distancia 1 que representa la distancia desde el centro de la tierra hasta el centro de la sonda"""


def CalculoDeDistanciaUno(x, y):
    d1 = math.sqrt((X1 - x) ** 2 + (Y1 - y) ** 2)
    return d1


def CalculoDeDistanciaDos(x, y):
    d2 = math.sqrt((X2 - x) ** 2 + (Y2 - y) ** 2)
    return d2


def CalculoDeDistanciaG(x, y):
    dg = math.sqrt(x ** 2 + y ** 2)
    return dg


def CalculoDeAlfaF1(x, y):
    alfaF1 = math.atan2(Y1 - y, X1 - x)
    return alfaF1


def CalculoDeAlfaF2(x, y):
    alfaF2 = math.atan2(Y2 - y, X2 - x)
    return alfaF2


def CalculoDeAlfaFc(x, y):
    alfaFc = math.atan2(y, x)
    return alfaFc


def SetearCondicionesIniciales():
    X.append(X1 - R1 - H0)
    Y.append(0)
    v0 = CalcularModuloVelocidadDeOrbitaCircular()
    alpha0 = math.pi / 2
    Vx.append(v0 * math.cos(alpha0))
    Vy.append(v0 * math.sin(alpha0))


def CalcularR():
    return R1 + H0


def CalcularModuloVelocidadDeOrbitaCircular():
    R = CalcularR()
    #v0 = math.sqrt((G * M1) / R)
    v0 = 10480.56
    return v0


def CalcularPeriodoT():
    R = CalcularR()
    v0 = CalcularModuloVelocidadDeOrbitaCircular()
    return (2 * math.pi * R) / v0


def CalculoDeDatos(x, y):
    d1 = CalculoDeDistanciaUno(x, y)
    d2 = CalculoDeDistanciaDos(x, y)
    dg = CalculoDeDistanciaG(x, y)
    alfa1 = CalculoDeAlfaF1(x, y)
    alfa2 = CalculoDeAlfaF2(x, y)
    alfac = CalculoDeAlfaFc(x, y)
    datos = (d1, d2, dg, alfa1, alfa2, alfac)
    return datos


def CalculoDeF_x(x, y):
    datos = CalculoDeDatos(x, y)
    f1 = ((G * M1) / datos[0] ** 2) * math.cos(datos[3])
    f2 = ((G * M2) / datos[1] ** 2) * math.cos(datos[4])
    fc = W ** 2 * datos[2] * math.cos(datos[5])
    return f1 + f2 + fc


def CalculoDeF_y(x, y):
    datos = CalculoDeDatos(x, y)
    f1 = ((G * M1) / datos[0] ** 2) * math.sin(datos[3])
    f2 = ((G * M2) / datos[1] ** 2) * math.sin(datos[4])
    fc = W ** 2 * datos[2] * math.sin(datos[5])
    return f1 + f2 + fc


def MetodoEuler():
    for i in range(N):
        Fxi = CalculoDeF_x(X[i], Y[i])
        Vx.append(Vx[i] + H * Fxi)

        Fyi = CalculoDeF_y(X[i], Y[i])
        Vy.append(Vy[i] + H * Fyi)

        X.append(X[i] + H * Vx[i])

        Y.append(Y[i] + H * Vy[i])


def MetodoRungeKutta2():
    for i in range(N):
        q1x = H * Vx[i]
        q1y = H * Vy[i]
        q1vx = H * CalculoDeF_x(X[i], Y[i])
        q1vy = H * CalculoDeF_y(X[i], Y[i])

        q2x = H * (Vx[i] + q1vx)
        q2y = H * (Vy[i] + q1vy)
        q2vx = H * CalculoDeF_x(X[i] + q1x, Y[i] + q1y)
        q2vy = H * CalculoDeF_y(X[i] + q1x, Y[i] + q1y)

        X.append(X[i] + (1 / 2) * (q1x + q2x))
        Y.append(Y[i] + (1 / 2) * (q1y + q2y))
        Vx.append(Vx[i] + (1 / 2) * (q1vx + q2vx))
        Vy.append(Vy[i] + (1 / 2) * (q1vy + q2vy))


SetearCondicionesIniciales()
MetodoRungeKutta2()

m = math.sqrt((X[0] - X[N - 1]) ** 2 + (Y[0] - Y[N - 1]) ** 2)
print("Diferencia Distancia:", m)


def CalculoEnergiaMecanica(i):
    v = math.sqrt(Vx[i] ** 2 + Vy[i] ** 2)
    d1 = CalculoDeDistanciaUno(X[i], Y[i])
    Ec = (1 / 2) * v ** 2
    Ep = -G * M1 / d1
    return Ec + Ep


print("Emecanica: ", CalculoEnergiaMecanica(N - 1) - CalculoEnergiaMecanica(0))

tierra = plt.Circle((X1, Y1), R1, color='blue')
luna = plt.Circle((X2, Y2), R2, color='grey')

fig, ax = plt.subplots()
plt.axis([-100000000, 100000000, -100000000, 100000000])
ax.add_artist(tierra)
ax.add_artist(luna)

plt.plot(X, Y)
plt.show()



