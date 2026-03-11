
# Problema 3: Estudio y visualización de la anomalía gravitatoria
# que genera una anomalía de masa en el subsuelo terrestre.


import numpy as np
import matplotlib.pyplot as plt


# Discretización de la superficie 5 km * 5 km

Dis = 200 # Divisiones realizadas para discretizar la superficie
dx = 5000 / Dis     #Longitud en x del vóxel
dy = 5000 / Dis
dz = 5000 / Dis
x = np.linspace(0, 5000, Dis)
y = np.linspace(0, 5000, Dis)
Xmalla, Ymalla = np.meshgrid(x, y)

d_prom = 1900 #Densidad promedio de suelo rocoso en kg/m^3

Discretizacion = {
    'Xmalla': Xmalla,
    'Ymalla': Ymalla,
    'dx': dx,
    'dy': dy,
    'dz': dz,
    'vol_celda' : dx * dy * dz,    # Volumen de cada celda en m^3
    'altura': 0,
    'd_suelo' : d_prom,
    'Dis' : Dis,

}


# DISCRETIZACIÓN DE LA ANOMALÍA DE MASA


# Anomalía de masa CÚBICA

config_cubo = {

    'd_escalar' : 2500,
    'longitud_lado' : 250,
    
    'centro' : (1500, 3000, 1000), #Posición del centro de la anomalía. (Consideramos el valor absoluto de z)
}

def crear_anomalia_cubo (config, superficie):

    dx = superficie ['dx']
    dy = superficie ['dy']
    dz = superficie ['dz']
    Dis = superficie ['Dis']
    d_cubo = np.zeros((Dis, Dis, Dis,))

    d_var = config ['d_escalar'] - superficie ['d_suelo']
    L = config ['longitud_lado']
    x0, y0, z0 = config ['centro']
    
    i_cx = int(x0 / dx)
    j_cy = int(y0 / dy)
    k_cz = int(z0 / dz)
    L_ind = int(L / dx)      #Traducimos los datos de metros a índices

    x_min = max(0, i_cx - L_ind // 2)
    x_max = min(Dis, i_cx + L_ind // 2)
    
    y_min = max(0, j_cy - L_ind // 2)
    y_max = min(Dis, j_cy + L_ind // 2)
    
    z_min = max(0, k_cz - L_ind // 2)
    z_max = min(Dis, k_cz + L_ind // 2)

    d_cubo [z_min : z_max, x_min : x_max, y_min : y_max] = d_var

    return {
        'densidad' : d_cubo,
        'vol_celda' : dx * dy * dz,
        'Dis' : Dis,
        'config' : config,
        'd_var' : d_var

    }




# Anomalía generada por cualquier PARALELEPÍPEDO


config_paralelepipedo = {

    'd_escalar' : 2800,
    'a' : 4000,    #Anchura, (sin girar se corresponde con x)
    'b' : 200,    #Altura, (sin girar se corresponde con y)
    'c' : 2000,     #Profundidad, (sin girar se corresponde con z)
    
    'centro' : (2500, 2500, 1500), #Posición del centro de la anomalía. (Consideramos el valor absoluto de z)
}


def crear_anomalia_paralelepipedo (config, superficie):

    dx = superficie ['dx']
    dy = superficie ['dy']
    dz = superficie ['dz']
    Dis = superficie ['Dis']
    d_paralel = np.zeros((Dis, Dis, Dis,))

    d_var = config ['d_escalar'] - superficie ['d_suelo']
    a, b, c = config ['a'], config ['b'], config ['c']
    x0, y0, z0 = config ['centro']
    
    i_cx = int(x0 / dx)
    j_cy = int(y0 / dy)
    k_cz = int(z0 / dz)
    a_ind = int(a / dx)
    b_ind = int(b / dy)
    c_ind = int(c / dz)     #Traducimos los datos de metros a índices

    x_min = max(0, i_cx -  a_ind// 2)
    x_max = min(Dis, i_cx + a_ind // 2)
    
    y_min = max(0, j_cy - b_ind // 2)
    y_max = min(Dis, j_cy + b_ind // 2)
    
    z_min = max(0, k_cz - c_ind // 2)
    z_max = min(Dis, k_cz + c_ind // 2)

    d_paralel [z_min : z_max, x_min : x_max, y_min : y_max] = d_var

    return {
        'densidad' : d_paralel,
        'vol_celda' : dx * dy * dz,
        'Dis' : Dis,
        'config' : config,
        'd_var' : d_var

    }



# Anomalía generada por una ESFERA

config_esfera = {

    'd_escalar' : 1000,
    'radio' : 500,
    
    'centro' : (2500, 2500, 1000),
}


def crear_anomalia_esfera (config, superficie):
    eje = np.linspace (0, 5000, superficie ['Dis'])
    Z, X, Y = np.meshgrid(eje, eje, eje, indexing ='ij')    #Introducimos la función indexing 'ij' porque meshgrid genera la malla en el orden X, Y, Z

    radio = config ['radio']   
    d_var = config ['d_escalar'] - superficie ['d_suelo']
    x0, y0, z0 = config ['centro']
    
    distancia_sq = (X - x0) ** 2 + (Y - y0) ** 2 + (Z - z0) ** 2
    es_legal = distancia_sq <= radio ** 2

    densidad_esfera = np.zeros ((superficie ['Dis'], superficie ['Dis'], superficie ['Dis']))
    densidad_esfera [es_legal] = d_var
    return {
        'densidad' : densidad_esfera,
        'vol_celda' : dx * dy * dz,
        'Dis' : Dis,
        'config' : config,
        'd_var' : d_var
                         
    }


#ROTACIÓN DE LA ANOMALÍA

phi = 0                 # Ángulo de giro alrededor del eje z
theta = 60               # Ángulo de giro alrededor del eje x
psi = 45                 # Ángulo de giro alrededor del nuevo eje z'



def matriz_rotacion (phi, theta, psi):
    rad_phi = np.radians(phi)
    rad_theta = np.radians(theta)
    rad_psi = np.radians(psi)
    D = np.array([
        [np.cos(rad_phi), np.sin(rad_phi), 0],
        [-np.sin(rad_phi), np.cos(rad_phi), 0],
        [0, 0, 1]    
    ])

    C = np.array([
        [1, 0, 0],
        [0, np.cos(rad_theta), -np.sin(rad_theta)],
        [0, -np.sin(rad_theta), np.cos(rad_theta)]
    ])

    B = np.array ([
        [np.cos(rad_psi), np.sin(rad_psi), 0],
        [-np.sin(rad_psi), np.cos(rad_psi), 0],
        [0, 0, 1]
    ])

    A = B @ C @ D
    return A


def rotacion_anomalia (anomalia, phi, theta, psi):
    R = matriz_rotacion(phi, theta, psi)
    eje = np.linspace (0, 5000, anomalia['Dis'])
    Z, X, Y = np.meshgrid(eje, eje, eje, indexing='ij')

    (x0, y0, z0) = anomalia ['config']['centro']

    X_rel = X - x0
    Y_rel = Y - y0
    Z_rel = Z - z0

    X_rot = X_rel * R[0,0] + Y_rel * R[0,1] + Z_rel * R[0,2]
    Y_rot = X_rel * R[1,0] + Y_rel * R[1,1] + Z_rel * R[1,2]
    Z_rot = X_rel * R[2,0] + Y_rel * R[2,1] + Z_rel * R[2,2]     #Sistema de coordenadas rotado

    radio = anomalia ['config']['longitud_lado'] / 2   #Aquí no hace falta transformar a índice de la matriz porque al ser meshgrid ya está en metros
    es_legal = (np.abs (X_rot) <= radio) & (np.abs (Y_rot) <= radio) & (np.abs (Z_rot) <= radio)  #Matriz (Dis, Dis, Dis) que contiene solo True y False
    densidad = np.zeros ( (anomalia ['Dis'], anomalia ['Dis'], anomalia ['Dis']) )
    densidad [es_legal] = anomalia ['d_var'] #En aquellos puntos donde es_legal tiene True, se asigna el valor de variación de densidad de la anomalía

    return densidad


# Rotación general para cualquier paralelepípedo

def rotacion_anomalia_general (anomalia, phi, theta, psi):
    R = matriz_rotacion(phi, theta, psi)
    eje = np.linspace (0, 5000, anomalia['Dis'])
    Z, X, Y = np.meshgrid(eje, eje, eje, indexing = 'ij')

    (x0, y0, z0) = anomalia ['config']['centro']

    X_rel = X - x0
    Y_rel = Y - y0
    Z_rel = Z - z0

    X_rot = X_rel * R[0,0] + Y_rel * R[0,1] + Z_rel * R[0,2]
    Y_rot = X_rel * R[1,0] + Y_rel * R[1,1] + Z_rel * R[1,2]
    Z_rot = X_rel * R[2,0] + Y_rel * R[2,1] + Z_rel * R[2,2]     #Sistema de coordenadas rotado

    radio_a = anomalia ['config']['a'] / 2
    radio_b = anomalia ['config']['b'] / 2
    radio_c = anomalia ['config']['c'] / 2
   #Aquí no hace falta transformar a índice de la matriz porque al ser meshgrid ya está en metros

    es_legal = (np.abs (X_rot) <= radio_a) & (np.abs (Y_rot) <= radio_b) & (np.abs (Z_rot) <= radio_c)  #Matriz (Dis, Dis, Dis) que contiene solo True y False
    densidad = np.zeros ( (anomalia ['Dis'], anomalia ['Dis'], anomalia ['Dis']) )
    densidad [es_legal] = anomalia ['d_var'] #En aquellos puntos donde es_legal tiene True, se asigna el valor de variación de densidad de la anomalía

    return densidad
         

#CÁLCULO DE LA ANOMALÍA GRAVITATORIA

def gravedad_vertical(anomalia, superficie):
    G = 6.67430e-11  # Constante gravitacional universal en m^3 kg^-1 s^-2
    indices_k, indices_i, indices_j = np.where (anomalia['densidad'] != 0)  #Coge solo los índices donde hay variación de densidad
    gz = np.zeros_like ((superficie ['Xmalla']))
  
    densidad = anomalia ['densidad'] [indices_k, indices_i, indices_j]
    masa = densidad * anomalia ['vol_celda']

    x_masa = indices_i * superficie ['dx']
    y_masa = indices_j * superficie ['dy']
    z_masa = indices_k * superficie ['dz']    #Transformación de los índices donde hay anomalía a metros reales
        
    for m, zv, xv, yv in zip(masa, z_masa, x_masa, y_masa):
        
        # Distancia desde ESTE bloque (xv, yv, zv) a TODA la superficie
        dist_X = superficie ['Xmalla'] - xv
        dist_Y = superficie['Ymalla'] - yv
        dist_Z = superficie ['altura'] + zv  #En python al hacer matriz - escalar, se resta a todas las componentes

        R = np.sqrt (dist_X**2 + dist_Y**2 + dist_Z**2)

        gz += (G * m * abs(dist_Z)) / (R**3)   # Se suma la gravedad en z por cada celda de la anomalía
        # Conversión a mGal (1m/s^2 = 100.000 mGal)
        gz_mGal = gz * 100000

    return gz_mGal  #Matriz como la superficie que en cada celda guarda la gravedad que aporta la anomalía

    


# RESULTADOS FINALES


esfera = crear_anomalia_esfera(config_esfera, Discretizacion)
cubo = crear_anomalia_cubo(config_cubo, Discretizacion)
dique = crear_anomalia_paralelepipedo(config_paralelepipedo, Discretizacion)

cubo ['densidad'] = rotacion_anomalia(cubo, phi, theta, psi)
dique ['densidad'] = rotacion_anomalia_general (dique, phi, theta, psi)


mapa_pot_cubo = gravedad_vertical(dique, Discretizacion)

L_total = 5000

plt.imshow(mapa_pot_cubo, extent = [0, L_total, 0, L_total], origin = 'lower', cmap = 'jet')
plt.colorbar(label='Anomalía Gravitatoria (mGal)')
plt.title('Mapa Gravimétrico generado por: cubo rotado')
plt.xlabel('Distancia X (m)')
plt.ylabel('Distancia Y (m)')
plt.grid(False) 

plt.show() 



