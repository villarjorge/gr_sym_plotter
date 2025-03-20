"""
Contiene las funciones para representar cosas en matplotlib. 
Queremos representar la superficie dada al tomar una variable de un cambio de coordenadas como 
constante (Ejemplo: en coordenadas esféricas tomamos constante el radio y tenemos una esfera)
"""
import matplotlib.pyplot as plt
import numpy as np

from symbolics import tomar_variables, get_lambdas_reducidas, get_lambda_escalar_curvatura

def plot_superficie(max_var1: float, max_var2: float, x_lambda, y_lambda, z_lambda):
    """
    Representa una superficie parametrizada con dos variables. La superficie tiene una serie de puntos (x, y, z) 
    que dependen de (u, v). x_lambda, y_lambda y z_lambda son el resultado de la función get_lamdas_reducidas

    Util: https://matplotlib.org/stable/gallery/mplot3d/surface3d_2.html#sphx-glr-gallery-mplot3d-surface3d-2-py
    https://matplotlib.org/stable/gallery/mplot3d/surface3d.html#sphx-glr-gallery-mplot3d-surface3d-py
    """
    u = np.linspace(0, max_var1, 100)
    v = np.linspace(0, max_var2, 100)

    U, V = np.meshgrid(u, v)

    X = x_lambda(U, V)
    Y = y_lambda(U, V)
    Z = z_lambda(U, V)

    # Plot the surface
    _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(X, Y, Z)

    # Set an equal aspect ratio
    ax.set_aspect('equal')

    plt.show()

def plot_superficie_color_r(max_var1: float, max_var2: float, escalar_lambda, x_lambda, y_lambda, z_lambda):
    """
    Representa una superficie parametrizada con dos variables. La superficie tiene una serie de puntos (x, y, z) 
    que dependen de (u, v). x_lambda, y_lambda y z_lambda son el resultado de la función get_lamdas_reducidas.

    La superficie viene coloreada por su escalar de curvatura

    https://likegeeks.com/coloring-python-3d-plots-matplotlib/
    """
    # Vectorizamos todas las lambdas porque es posible que si no de errores
    vec_escalar_lambda, vec_x_lambda, vec_y_lambda, vec_z_lambda = np.vectorize(escalar_lambda), np.vectorize(x_lambda), np.vectorize(y_lambda), np.vectorize(z_lambda)

    u = np.linspace(0, max_var1, 100)
    v = np.linspace(0, max_var2, 100)

    U, V = np.meshgrid(u, v)

    X = vec_x_lambda(U, V)
    Y = vec_y_lambda(U, V)
    Z = vec_z_lambda(U, V)

    E = vec_escalar_lambda(U, V)

    # Plot the surface
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    superficie = ax.plot_surface(X, Y, Z, facecolors=plt.cm.viridis(E)) # np.nan_to_num(E) Reemplazar nans con valores grandes

    # Añadimos una colorbar
    fig.colorbar(superficie)

    # Set an equal aspect ratio
    ax.set_aspect('equal')

    plt.show()


if __name__ == "__main__":
    print("Testando representaciones gráficas")

    # variables_todas = tomar_variables("r theta phi")

    # lambdas = get_lambdas_reducidas(variables_todas, "r*sin(theta)*cos(phi)", "r*sin(theta)*sin(phi)", "r*cos(theta)", 0, 1) # R = 1 (coord 0 igual a uno)
    # print("Testando representación sin coloreado")
    # plot_superficie(np.pi, 2*np.pi, *lambdas)

    # print("Testando representación con coloreado")
    # escalar_lambda = get_lambda_escalar_curvatura(variables_todas, "r*sin(theta)*cos(phi)", "r*sin(theta)*sin(phi)", "r*cos(theta)", 0, 1)
    # plot_superficie_color_r(np.pi, 2*np.pi, escalar_lambda, *lambdas)

    print("Probando con toroide")

    b = 2
    x_str, y_str, z_str = f"({b} + r*cos(phi))*cos(theta)", f"({b} + r*cos(phi))*sin(theta)", f"r*sin(phi)"

    variables_todas = tomar_variables("r theta phi")

    coord, valor = 0, 1

    lambdas = get_lambdas_reducidas(variables_todas, x_str, y_str, z_str, coord, valor) # R = 1 (coord 0 igual a uno)
    escalar_lambda = get_lambda_escalar_curvatura(variables_todas, x_str, y_str, z_str, coord, valor)

    # print("Testando representación sin coloreado")
    # plot_superficie(2*np.pi, 2*np.pi, *lambdas)

    print("Testando representación con coloreado")
    plot_superficie_color_r(2*np.pi, 2*np.pi, escalar_lambda, *lambdas)

    # Hiperboloide

    # x_str, y_str, z_str = "r*cosh(phi)*cos(theta)", "r*cosh(phi)*sin(theta)", "r*sinh(phi)"

    # variables_todas = tomar_variables("r theta phi")

    # coord, valor = 0, 1

    # lambdas = get_lambdas_reducidas(variables_todas, x_str, y_str, z_str, coord, valor) # R = 1 (coord 0 igual a uno)
    # escalar_lambda = get_lambda_escalar_curvatura(variables_todas, x_str, y_str, z_str, coord, valor)

    # print("Testando representación con coloreado")
    # plot_superficie_color_r(10, 2*np.pi, escalar_lambda, *lambdas)