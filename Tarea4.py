#######################################################
'''
Metodos Numericos para la Ciencia e Ingenieria
FI3104-1
Tarea 3
Maximiliano Dirk Vega Aguilera
18.451.231-9
'''
#######################################################

import numpy as np
import matplotlib.pyplot as plt
from planeta import Planeta

#######################################################

def ec_mov(p, t, dt):
    '''
    Integra la ecuacion de movimiento utilizando euler, rk4 o verlet
    '''
    x = np.zeros(len(t))
    y = np.zeros(len(t))
    E = np.zeros(len(t))
    for i in range(len(t)):
        x[i] = p.y_actual[0]
        y[i] = p.y_actual[1]
        t[i] = p.t_actual
        E[i] = p.energia_total()
        #p.avanza_euler(dt)   #cambiar p.avanza_* por metodo deseado
        p.avanza_rk4(dt)
        #p.avanza_verlet(dt)

    return x, y, t, E

def perihelio(x, y):
    '''
    determina la posicion del perihelio mediante la comparacion del radio en
    los distintos pasos
    '''
    peri_x = np.array([])
    peri_y = np.array([])
    for i in range(len(x)-2):
        if x[i+1]**2 + y[i+1]**2 >= x[i]**2 + y[i]**2:
            if x[i+1]**2 + y[i+1]**2 >= x[i+2]**2 + y[i+2]**2:
                peri_x = np.append(peri_x, x[i+1])
                peri_y = np.append(peri_y, y[i+1])

    return peri_x, peri_y

#condiciones iniciales
vyo = 0.4  #0.4 para euler y rk4  0.05 para verlet
condicion_inicial = [10., 0., 0., vyo]
alpha = 10**(-2.231)  #231 ultimos 3 digitos antes de verificador

p = Planeta(condicion_inicial,alpha)  #agregar o quitar alpha para los casos 2 y 3

#se definen variables a usar
n = 1000000    #numero de pasos a realizar  #10000 sin alpha....100000 con alpha verlet
                                                            #   1000000 con alpha euler y rk4
t  = np.zeros(n)  #genera vector tiempo a completar
dt = 0.05         #tamanho del paso
x  = np.zeros(n)  #genera vector posicion x a completar
y  = np.zeros(n)  #genera vector posicion y a completar
E = np.zeros(n)   #genera vector energia total a completar

#se completan las variables para graficar
x, y, t, E = ec_mov(p, t, dt)
peri_x, peri_y = perihelio(x,y)

print 'n orbitas =', len(peri_x) #identifica numero de orbitaciones

plt.plot(x,y)
#plt.title('Orbita con Euler')
plt.title('Orbita con Runge-Kutta 4')
#plt.title('Orbita con Verlet')
plt.xlabel("x")
plt.ylabel("y")
plt.show()

plt.plot(t, E)
#plt.title('Energia total del sistema con Euler')
plt.title('Energia total del sistema con Runge-Kutta 4')
#plt.title('Energia total del sistema con Verlet')
plt.xlabel("tiempo")
plt.ylabel("energia")
plt.show()

plt.plot(x,y,peri_x,peri_y,'o')
#plt.title('Ubicacion del perihelio en la Orbita con Euler')
plt.title('Ubicacion del perihelio en la Orbita con Runge-Kutta 4')
#plt.title('Ubicacion del perihelio en la Orbita con con Verlet')
plt.xlabel("x, perihelio")
plt.ylabel("y")
plt.show()
