"""
En este archivo se encuentran las funciones que calculan simbólicamente la métrica, 
simbolos de Cristoefel etc, así como funciones helper que facilitan el paso de variables 
entre las funciones.

Notas para el futuro: 
    - tal vez se podrían combinar las funciones get_lambdas_reducidas y get_lambda_escalar_curvatura
    - En general me gustaría que las funciones fuesen más generales (infiriesen variables?), mira testing2
    - ¿Son necesareas las variables que no dependen de tau?
"""
import sympy as smp
from itertools import product

def tomar_variables(variables: str):
    """
    Toma las nuevas variables como string y devuelve una lista de simbolos de sympy
    Crea tanto variables que dependen de tau como que no dependen de tau, y tau
    """
    tau = smp.symbols("tau")

    u, v, w = smp.symbols(variables, cls=smp.Function)
    
    u = u(tau)
    v = v(tau)
    w = w(tau)

    return (u, v, w), smp.symbols(variables)

def get_lambdas(variables_todas, x: str, y: str, z: str):
    """
    Convierte el cambio de coordenadas simbólico a uno numérico
    """
    _, variables = variables_todas

    x, y, z = smp.sympify(x), smp.sympify(y), smp.sympify(z)

    x_lambda = smp.lambdify(variables, x)
    y_lambda = smp.lambdify(variables, y)
    z_lambda = smp.lambdify(variables, z)

    return x_lambda, y_lambda, z_lambda

def get_lambdas_reducidas(variables_todas, x: str, y: str, z: str, m: int, valor: float):
    """
    Versión reducida de las lambdas donde la coordenada m se ha sustituido por el valor "valor"
    """
    _, variables = variables_todas

    x, y, z = smp.sympify(x), smp.sympify(y), smp.sympify(z)
    x, y, z = x.subs(variables[m], valor), y.subs(variables[m], valor), z.subs(variables[m], valor)

    variables = list(variables)
    variables.pop(m)

    x_lambda = smp.lambdify(variables, x)
    y_lambda = smp.lambdify(variables, y)
    z_lambda = smp.lambdify(variables, z)

    return x_lambda, y_lambda, z_lambda

def calcular_metrica(variables_todas, x: str, y: str, z: str):
    """
    Calcula la métrica dado un cambio de coordenadas x, y, z = f(variables), g(variables), h(variables)
    variables debe ser algo como "u v w" o r"r theta phi" 

    Ejemplo: Métrica del cambio a coords esféricas. variables = lista se variables de sympy, usar tomar variables;
    x = "r*sin(theta)*cos(phi)"; y = "r*sin(theta)*sin(phi)"; z = "r*cos(theta)"

    Util: https://physics.stackexchange.com/questions/321781/how-is-the-spherical-coordinate-metric-tensor-derived
    """
    variables_tau, variables = variables_todas # Primeras dependen de tau y segundas no

    x, y, z = smp.sympify(x), smp.sympify(y), smp.sympify(z)

    g = [[0]*3 for _ in range(3)] # this does not work: [[0]*3]*3

    for i in range(0, 3):
        for j in range(0, 3):
            for v in [x, y, z]:
                # Si derivasemos v por theta(tau) nos daría cero, ya que para el programa theta no es theta(tau)
                v = v.subs(variables[0], variables_tau[0]).subs(variables[1], variables_tau[1]).subs(variables[2], variables_tau[2])
                # Implicitamente g_(ij) es la delta de kroneker
                g[i][j] += smp.diff(v, variables_tau[i])*smp.diff(v, variables_tau[j])
            g[i][j] = g[i][j].simplify()

    return g


def reducir_metrica(g, variables_todas, m: int):
    """
    Quita la dimensión m de la métrica g dada, lo que reduce las coordenadas 
    a una superficie (esto es un embebimiento)
    variables es lo que retorna tomar_variables
    m debe ser 0, 1 ó 2
    """
    variables_tau, variables = variables_todas # Primeras dependen de tau y segundas no

    new_g = [[0]*2 for _ in range(2)]

    indexer = [0, 1, 2]
    indexer.remove(m)

    for i in range(0, 2):
        for j in range(0, 2):
            new_g[i][j] = g[indexer[i]][indexer[j]].subs(variables_tau[m], variables[m])

    return new_g

def imprimir_metrica(g):
    for i in range(0, len(g)):
        for j in range(0, len(g)):
            print(g[i][j], end="  ")
        print()

def calcular_eel(g, variables_todas, m: int):
    """
    Calcula el lagrangiano de la métrica y calcula las ecuaciones de euler lagrange (eel)
    variables es lo que retorna tomar_variables
    """
    variables_tau, _ = variables_todas # Primeras dependen de tau y segundas no

    indexer = [0, 1, 2]
    indexer.remove(m)

    tau = smp.symbols("tau")

    L = smp.sympify(0)

    for i in range(0, 2):
        for j in range(0, 2):
            L += g[i][j]*smp.diff(variables_tau[indexer[i]], tau)*smp.diff(variables_tau[indexer[j]], tau)

    # Estas tres lineas toman algo de tiempo
    e1, e2 = smp.euler_equations(L, variables_tau, tau)
    # e1 = e1.simplify() # Comentadas por ahora 
    # e2 = e2.simplify()

    return e1, e2

def calcular_christoffel(eel, variables_todas, m: int):
    Gamma = [[[0 for _ in range(2)] for _ in range(2)] for _ in range(2)] # Matriz 2*2*2

    variables_tau, _ = variables_todas # Primeras dependen de tau y segundas no

    indexer = [0, 1, 2]
    indexer.remove(m)

    tau = smp.symbols("tau")

    for i in range(0, 2): # primer indice
        # Resolvemos d^2/dt^2 (v^i) = elem
        elem = smp.solve(eel[i], smp.diff(variables_tau[indexer[i]], tau, tau))[0]
        elem = smp.expand(elem)
        # elem será una suma de coeficientes lineares en las derivadas de las variables_tau
        # Por ahora hacemos un loop doble, pero tal vez sería posible extraer las derivadas directamente
        # https://docs.sympy.org/latest/index.html
        for j in range(0, 2):
            for k in range(0, 2):
                var = smp.diff(variables_tau[indexer[j]], tau)*smp.diff(variables_tau[indexer[k]], tau)
                Gamma[i][j][k] = -elem.coeff(var)
                if j != k: # The common pitfall
                    Gamma[i][j][k] = Gamma[i][j][k]/2

    return Gamma

def imprimir_christoffel(Gamma, variables_todas, m: int):
    """
    Imprime los símbolos de christoffel no ceros
    """
    _, vari = variables_todas # Primeras dependen de tau y segundas no

    indexer = [0, 1, 2]
    indexer.remove(m)
    
    for i, j, k in product(range(0, 2), repeat=3):
        if Gamma[i][k][j] != 0:
            print(f"Gamma^{vari[indexer[i]]}_({vari[indexer[j]]} {vari[indexer[k]]}) = {Gamma[i][j][k]}")

def calcular_riemman(Gamma, variables_todas, m: int):
    """
    Calcula el tensor de riemman de manera directa 
    """
    Riemman = [[[[0 for _ in range(2)] for _ in range(2)] for _ in range(2)] for _ in range(2)] # Matriz 2*2*2*2

    variables_tau, _ = variables_todas # Primeras dependen de tau y segundas no

    indexer = [0, 1, 2]
    indexer.remove(m)
    
    for i, j, k, l in product(range(0, 2), repeat=4):
        Riemman[i][j][k][l] = smp.diff(Gamma[i][l][j], variables_tau[indexer[k]]) - smp.diff(Gamma[i][k][j], variables_tau[indexer[l]])
        for alpha in range(0, 2):
            Riemman[i][j][k][l] += Gamma[i][k][alpha]*Gamma[alpha][l][j] - Gamma[i][l][alpha]*Gamma[alpha][k][j]
        # Simplificaciones para la métrica de la esfera, simplificar otras métricas requerirá en general otras operaciones
        # https://docs.sympy.org/latest/tutorials/intro-tutorial/simplification.html#trigsimp
        Riemman[i][j][k][l] = smp.expand_trig(Riemman[i][j][k][l])
        Riemman[i][j][k][l] = smp.trigsimp(Riemman[i][j][k][l])
        # Riemman[i][j][k][l] = smp.expand(Riemman[i][j][k][l])
        # Riemman[i][j][k][l] = smp.simplify(Riemman[i][j][k][l])

    return Riemman

def imprimir_riemman(Riemman, variables_todas, m: int):
    """
    Imprime los símbolos de christoffel no ceros
    """
    _, vari = variables_todas # Primeras dependen de tau y segundas no

    indexer = [0, 1, 2]
    indexer.remove(m)
    
    for i, j, k, l in product(range(0, 2), repeat=4):
        if Riemman[i][j][k][l] != 0:
            print(f"Riemman^{vari[indexer[i]]}_({vari[indexer[j]]} {vari[indexer[k]]} {vari[indexer[l]]}) = {Riemman[i][j][k][l]}")

def calcular_ricci(Riemman):
    """
    Contrae el tensor de Riemman para encontrar el tensor de Ricci
    """
    Ricci = [[0]*2 for _ in range(2)] # this does not work: [[0]*3]*3

    for i, j in product(range(0, 2), repeat=2):
        for alpha in range(0, 2):
            Ricci[i][j] += Riemman[alpha][i][alpha][j]

    return Ricci

def imprimir_ricci(Ricci, variables_todas, m):
    _, vari = variables_todas # Primeras dependen de tau y segundas no

    indexer = [0, 1, 2]
    indexer.remove(m)

    for i, j in product(range(0, 2), repeat=2):
        if Ricci[i][j] != 0:
            print(f"Ricci^{vari[indexer[i]]}_{vari[indexer[j]]} = {Ricci[i][j]}")

def calcular_escalar_curvatura(Ricci, g):
    """
    Contrae el tensor de Ricci para encontrar el escalar de curvatura.
    """
    R = smp.sympify(0)

    for i, j in product(range(0, 2), repeat=2):
        R += (smp.Matrix(g)**(-1)).row(i)[j]*Ricci[i][j]
    
    return R

def get_lambda_escalar_curvatura(variables_todas, x: str, y: str, z: str, m: int, valor: float):
    """
    Como get lambdas_reducidas, obtén una lambda del escalar de curvatura, con la coordenada m igual a "valor"
    """
    g_esfericas = calcular_metrica(variables_todas, x, y, z)
    g_red = reducir_metrica(g_esfericas, variables_todas, m)
    ec = calcular_eel(g_red, variables_todas, m)
    G = calcular_christoffel(ec, variables_todas, m)
    R = calcular_riemman(G, variables_todas, m)
    Ric = calcular_ricci(R)
    escalar = calcular_escalar_curvatura(Ric, g_red)

    variables_tau, variables = variables_todas

    # Reemplaza las variables en tau por variables que no dependen de tau

    for alpha in range(0, 3):
        escalar = escalar.subs(variables_tau[alpha], variables[alpha])

    escalar = escalar.subs(variables[m], valor)

    variables = list(variables)
    variables.pop(m)

    return smp.lambdify(variables, escalar)

if __name__ == "__main__":
    print("Probando que funcione, calculando la métrica en coordenadas esfericas")

    variables_sym = tomar_variables("r theta phi")
    print(variables_sym)
    # esfera
    # x_str, y_str, z_str = "r*sin(theta)*cos(phi)", "r*sin(theta)*sin(phi)", "r*cos(theta)"
    # Toro
    b = 2
    x_str, y_str, z_str = f"({b} + r*cos(phi))*cos(theta)", f"({b} + r*cos(phi))*sin(theta)", f"r*sin(phi)"

    g_esfericas = calcular_metrica(variables_sym, x_str, y_str, z_str)
    
    imprimir_metrica(g_esfericas)
    
    print()
    print("Probando a reducir la métrica")
    print()

    m = 0

    g_red = reducir_metrica(g_esfericas, variables_sym, m)

    imprimir_metrica(g_red)

    print()
    print("Calculando las EEL")
    print()

    ec = calcular_eel(g_red, variables_sym, m)

    print(ec[0])
    print(ec[1])

    # print(ec[0].args[0].as_coeff_add()[1])
    # print(ec[1].args[0])

    # tau = smp.symbols("tau")

    # test = smp.solve(ec[0], smp.diff(variables_sym[0][1], tau, tau)) # Resolver de esta manera compacta sin(theta)cos(theta)
    # print(test)
    # print(test.as_coefficient(smp.diff(variables_sym[0][1], tau)))
    # print(test.coeff(smp.diff(variables_sym[0][1], tau)))

    print()

    G = calcular_christoffel(ec, variables_sym, m)

    print("Imprimiendo simbolos de christoffel no nulas")
    imprimir_christoffel(G, variables_sym, m)
    print()

    R = calcular_riemman(G, variables_sym, m)
    print("Imprimiendo componentes del tensor de Riemman no nulas")
    imprimir_riemman(R, variables_sym, m)
    print()

    Ri = calcular_ricci(R)
    print("Imprimiendo componentes del tensor de Ricci no nulas")
    imprimir_ricci(Ri, variables_sym, m)
    print()

    print(f"Escalar de curvatura: {calcular_escalar_curvatura(Ri, g_red)}")