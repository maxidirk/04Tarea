#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

#######################################################

class Planeta(object):
    '''
    Complete el docstring.
    '''

    def __init__(self, condicion_inicial, alpha=0):
        '''
        __init__ es un método especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_anterior = condicion_inicial
        self.y_actual = condicion_inicial
        self.t_anterior = 0.
        self.t_actual = 0.
        self.alpha = alpha

    def ecuacion_de_movimiento(self,a=[0,0,0,0]): #arreglo auxiliar por problemas en rk4
        '''
        Implementa la ecuación de movimiento, como sistema de ecuaciónes de
        primer orden.
        Se considera G=M=m=1
        '''
        if a == [0,0,0,0]:
            a = self.y_actual
        x, y, vx, vy = a
        fx = (-x) * (1./np.sqrt(x**2 + y**2) - 2. * self.alpha/(x**2 + y**2)) / (x**2 + y**2)
        fy = (-y) * (1./np.sqrt(x**2 + y**2) - 2. * self.alpha/(x**2 + y**2)) / (x**2 + y**2)

        return [vx, vy, fx, fy]

    def avanza_euler(self, dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.
        '''
        x1, y1, vx1, vy1 = self.y_actual
        vx, vy, fx, fy = self.ecuacion_de_movimiento()
        x2  = x1  + vx * dt
        y2  = y1  + vy * dt
        vx2 = vx1 + fx * dt
        vy2 = vy1 + fy * dt
        self.y_anterior = self.y_actual
        self.y_actual   = [x2, y2, vx2, vy2]
        self.t_anterior = self.t_actual
        self.t_actual += dt

        pass

    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        '''
        x1, y1, vx1, vy1 = self.y_actual
        vx, vy, fx, fy = self.ecuacion_de_movimiento()
        k1 = dt * np.asarray(self.ecuacion_de_movimiento([x1,y1,vx,vy]))
        k2 = dt * np.asarray(self.ecuacion_de_movimiento([x1 + k1[0]/2., y1 + k1[1]/2. , vx + k1[2]/2., vy + k1[3]/2.]))
        k3 = dt * np.asarray(self.ecuacion_de_movimiento([x1 + k2[0]/2., y1 + k2[1]/2. , vx + k2[2]/2., vy + k2[3]/2.]))
        k4 = dt * np.asarray(self.ecuacion_de_movimiento([x1 + k3[0], y1 + k3[1] , vx + k3[2], vy + k3[3]]))
        x2  = x1 + k1[0]/6. + k2[0]/3. + k3[0]/3. + k4[0]/6.
        y2  = y1 + k1[1]/6. + k2[1]/3. + k3[1]/3. + k4[1]/6.
        vx2 = vx + k1[2]/6. + k2[2]/3. + k3[2]/3. + k4[2]/6.
        vy2 = vy + k1[3]/6. + k2[3]/3. + k3[3]/3. + k4[3]/6.
        self.y_anterior = self.y_actual
        self.y_actual = [x2, y2, vx2, vy2]
        self.t_anterior = self.t_actual
        self.t_actual += dt

        pass

    def avanza_verlet(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet.
        '''
        x0, y0, vx0, vy0 = self.y_anterior
        x1, y1, vx1, vy1 = self.y_actual
        vx, vy, fx, fy   = self.ecuacion_de_movimiento()
        if self.t_actual == 0:  #el primer valor se obtiene mediante rk4
            self.avanza_rk4(dt)
        else:                   #ya teniendo 2 valores anteriores, se procede con Verlet
            x2  = 2. * x1 - x0 + (0.5 * fx * dt**2)
            y2  = 2. * y1 - y0 + (0.5 * fy * dt**2)
            vx2 = (x2 - x1) / dt
            vy2 = (y2 - y1) / dt
            self.y_anterior = self.y_actual
            self.y_actual   = [x2, y2, vx2, vy2]
            self.t_anterior = self.t_actual
            self.t_actual += dt

        pass

    def energia_total(self):
        '''
        Calcula la enérgía total del sistema en las condiciones actuales.
        G=M=m=1
        '''
        x, y, vx, vy = self.y_actual
        T = 0.5 * (vx**2 + vy**2)  #energia cinetica
        V = (-1.) / np.sqrt(x**2 + y**2) + self.alpha / (x**2 + y**2) #energia potencial
        E = T + V #energia total

        return E

#######################################################
