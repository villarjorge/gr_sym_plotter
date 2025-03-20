"""

Nota: creo que la manera en la que tomo valor absoluto está mal, mirar diff_eqs2.py

Lo mismo que diffrential_eqs.ipynb pero en un archivo de python para ver mejor la gráfica
"""

import numpy as np
from scipy.integrate import odeint # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
import matplotlib.pyplot as plt

def toro(Y, t, a:float, b:float, L:float, l:float):
    phi, theta = Y
    w = (b+a*np.sin(phi))
    # Debemos devolver la derivada de Y, dY/dt = (dotphi, dottheta)
    # return [(1/a)*np.sqrt(L*w**2 - l**2)/w, l/w**2]
    # Valor absoluto ad hoc, ¿Es como añadir una const. al lagrangiano para que simepre sea positivo?
    return [(1/a)*np.sqrt(abs(L*w**2 - l**2))/w, l/w**2]

def resolver(t, a:float, b:float, phi_0:float, dotphi_0:float, theta_0:float, dottheta_0:float):
    # t es un linspace de los tiempos que quieras
    w = (b+a*np.sin(phi_0))
    l = dottheta_0*w**2
    L = (w**2)*(dottheta_0**2) + (a**2)*(dotphi_0)

    return odeint(toro, [phi_0, theta_0], t, args=(a, b, L, l))

def representar_geodesica(t, a, b, phi_0, dotphi_0, theta_0, dottheta_0):
    u = np.linspace(0, 2*np.pi, 100) # phi
    v = np.linspace(0, 2*np.pi, 100) # theta

    U, V = np.meshgrid(u, v)

    X = (b+a*np.sin(U))*np.cos(V)
    Y = (b+a*np.sin(U))*np.sin(V)
    Z = a*np.cos(U)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 8))

    superficie = ax.plot_surface(X, Y, Z, alpha=0.5)
    #superficie = ax.plot_wireframe(X, Y, Z)

    # Set an equal aspect ratio
    ax.set_aspect('equal')

    sol = resolver(t, a, b, phi_0=phi_0, dotphi_0=dotphi_0, theta_0=theta_0, dottheta_0=dottheta_0)
    phi = sol[:, 0]
    theta = sol[:, 1]

    x = (b+a*np.sin(phi))*np.cos(theta)
    y = (b+a*np.sin(phi))*np.sin(theta)
    z = a*np.cos(phi)

    geodesica = ax.plot(x, y, z, color="red", linewidth=3)

    x_0 = (b+a*np.sin(phi_0))*np.cos(theta_0)
    y_0 = (b+a*np.sin(phi_0))*np.sin(theta_0)
    z_0 = a*np.cos(phi_0)

    punto_inicial = ax.scatter(x_0, y_0, z_0, color="green")

    plt.title(rf"{a = }, {b = }, $\phi_0$ = {round(phi_0, 3)}, $\dot\phi_0$ = {round(dotphi_0, 3)}, $\theta_0$ = {round(theta_0, 3)}, $\dot\theta_0$ = {round(dottheta_0, 3)}")

    plt.show()

if __name__ == "__main__":
    # Geodesicas con ángulos constantes
    # La de phi constante
    # Sabemos que phi = 0, pi, pero aquí solo salen las que deberían cuando phi = pi/2, 3*pi/2
    representar_geodesica(t=np.linspace(0, 2*np.pi, 100), a=0.5, b=1, phi_0=3*np.pi/2, dotphi_0=0, theta_0=0, dottheta_0=1)
    # Las de theta constante
    # representar_geodesica(t=np.linspace(0, 2*np.pi, 100), a=0.5, b=1, phi_0=0, dotphi_0=1, theta_0=0, dottheta_0=0)

    # Otras geodesicas
    # representar_geodesica(t=np.linspace(0, 10*np.pi, 500), a=0.5, b=1, phi_0=np.pi/2, dotphi_0=2**(-1/2), theta_0=0, dottheta_0=2**(-1/2))
    # representar_geodesica(t=np.linspace(0, 10*np.pi, 1000), a=0.5, b=1, phi_0=np.pi/2, dotphi_0=1/2, theta_0=0, dottheta_0=(3/4)**0.5)