"""
Estas sé que están bien, en el otro pasan cosas raras con la raiz cuadrada
"""

import numpy as np
from scipy.integrate import odeint # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

def toro(Y, t, a:float, b:float):
    phi, dotphi, theta, dottheta = Y
    w = (b+a*np.sin(phi))
    # Devolvemos la derivada de Y: dY/dt
    return [dotphi, (np.cos(phi)*w/a)*dottheta**2, dottheta, (-2*a*np.cos(phi)/w)*dotphi*dottheta]

def resolver(t, a:float, b:float, phi_0:float, dotphi_0:float, theta_0:float, dottheta_0:float):
    # t es un linspace de los tiempos que quieras
    return odeint(toro, [phi_0, dotphi_0, theta_0, dottheta_0], t, args=(a, b))

def representar_geodesica(t, a, b, phi_0, dotphi_0, theta_0, dottheta_0):
    u = np.linspace(0, 2*np.pi, 100) # phi
    v = np.linspace(0, 2*np.pi, 100) # theta

    U, V = np.meshgrid(u, v)

    X = (b+a*np.sin(U))*np.cos(V)
    Y = (b+a*np.sin(U))*np.sin(V)
    Z = a*np.cos(U)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 8))

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    superficie = ax.plot_surface(X, Y, Z, alpha=0.5)
    #superficie = ax.plot_wireframe(X, Y, Z)

    # Set an equal aspect ratio
    ax.set_aspect('equal')

    sol = resolver(t, a, b, phi_0=phi_0, dotphi_0=dotphi_0, theta_0=theta_0, dottheta_0=dottheta_0)
    phi = sol[:, 0]
    # dotphi = sol[:, 1]
    theta = sol[:, 2]
    # dottheta = sol[:, 3]

    x = (b+a*np.sin(phi))*np.cos(theta)
    y = (b+a*np.sin(phi))*np.sin(theta)
    z = a*np.cos(phi)

    geodesica = ax.plot(x, y, z, color="red", linewidth=1)

    x_0 = (b+a*np.sin(phi_0))*np.cos(theta_0)
    y_0 = (b+a*np.sin(phi_0))*np.sin(theta_0)
    z_0 = a*np.cos(phi_0)

    punto_inicial = ax.scatter(x_0, y_0, z_0, color="green")

    plt.title(rf"{a = }, {b = }, $\phi_0$ = {round(phi_0, 3)}, $\dot\phi_0$ = {round(dotphi_0, 3)}, $\theta_0$ = {round(theta_0, 3)}, $\dot\theta_0$ = {round(dottheta_0, 3)}")
    plt.show()


def representar_geodesica_sliders(t, a, b):
    phi_0, dotphi_0, theta_0, dottheta_0 = 0, 1, 0, 0

    u = np.linspace(0, 2*np.pi, 100) # phi
    v = np.linspace(0, 2*np.pi, 100) # theta

    U, V = np.meshgrid(u, v)

    X = (b+a*np.sin(U))*np.cos(V)
    Y = (b+a*np.sin(U))*np.sin(V)
    Z = a*np.cos(U)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 8))

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    superficie = ax.plot_surface(X, Y, Z, alpha=0.5)
    #superficie = ax.plot_wireframe(X, Y, Z)

    # Set an equal aspect ratio
    ax.set_aspect('equal')

    axcolor = 'lightgoldenrodyellow' # definimos el color de los sliders que usaremos más adelante
    ax.margins(x=0)

    sol = resolver(t, a, b, phi_0=phi_0, dotphi_0=dotphi_0, theta_0=theta_0, dottheta_0=dottheta_0)
    phi = sol[:, 0]
    theta = sol[:, 2]

    x = (b+a*np.sin(phi))*np.cos(theta)
    y = (b+a*np.sin(phi))*np.sin(theta)
    z = a*np.cos(phi)

    geodesica = ax.plot(x, y, z, color="red", linewidth=1)

    x_0 = (b+a*np.sin(phi_0))*np.cos(theta_0)
    y_0 = (b+a*np.sin(phi_0))*np.sin(theta_0)
    z_0 = a*np.cos(phi_0)

    punto_inicial = ax.scatter(x_0, y_0, z_0, color="green")

    # ajustamos la posición de la representación para acomodar los sliders
    plt.subplots_adjust(left=0.30, bottom=0.30)

    # hacer un slider vertical para controlar "phi_0". Los cuatro números controlan la posición, el grosor y la altura del slider
    axamp = plt.axes([0.1, 0.25, 0.0225, 0.63], facecolor=axcolor)
    phi_0_slider = Slider(
        ax=axamp,
        label=r"$\phi_0$",
        valmin=0,
        valmax=np.pi*2,
        valinit=phi_0,
        orientation="vertical"
    )

    axamp = plt.axes([0.15, 0.25, 0.0225, 0.63], facecolor=axcolor)
    theta_0_slider = Slider(
        ax=axamp,
        label=r"$\theta_0$",
        valmin=0,
        valmax=np.pi*2,
        valinit=theta_0,
        orientation="vertical"
    )

    axamp = plt.axes([0.20, 0.25, 0.0225, 0.63], facecolor=axcolor)
    dotphi_0_slider = Slider(
        ax=axamp,
        label=r"$\dot\phi_0$",
        valmin=0,
        valmax=1,
        valinit=dotphi_0,
        orientation="vertical"
    )
    axamp = plt.axes([0.25, 0.25, 0.0225, 0.63], facecolor=axcolor)
    dottheta_0_slider = Slider(
        ax=axamp,
        label=r"$\dot\theta_0$",
        valmin=0,
        valmax=1,
        valinit=dottheta_0,
        orientation="vertical"
    )
    # La función que será llamada cuando se actualice el slider
    def update(val): # para cada vez que se mueve el slider recalculamos la solución
        ax.clear()

        superficie = ax.plot_surface(X, Y, Z, alpha=0.5)

        phi_0, dotphi_0, theta_0, dottheta_0 = phi_0_slider.val, dotphi_0_slider.val, theta_0_slider.val, dottheta_0_slider.val
        sol = resolver(t, a, b, phi_0=phi_0, dotphi_0=dotphi_0, theta_0=theta_0, dottheta_0=dottheta_0)
        phi = sol[:, 0]
        theta = sol[:, 2]

        x = (b+a*np.sin(phi))*np.cos(theta)
        y = (b+a*np.sin(phi))*np.sin(theta)
        z = a*np.cos(phi)

        geodesica = ax.plot(x, y, z, color="red", linewidth=1)

        x_0 = (b+a*np.sin(phi_0))*np.cos(theta_0)
        y_0 = (b+a*np.sin(phi_0))*np.sin(theta_0)
        z_0 = a*np.cos(phi_0)

        punto_inicial = ax.scatter(x_0, y_0, z_0, color="green")
    
    # Estas lineas registran cuando se cambia el slider
    phi_0_slider.on_changed(update)
    dotphi_0_slider.on_changed(update)
    theta_0_slider.on_changed(update)
    dottheta_0_slider.on_changed(update)
    
    # Creamos un boton para resetear los valores del slider
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    def reset(event): # esta función se llamará cuando se pulse el botón
        phi_0_slider.reset()
        dotphi_0_slider.reset()
        theta_0_slider.reset()
        dottheta_0_slider.reset()
    button.on_clicked(reset)

    plt.show()


def representar_geodesica_sliders2(t, a, b):
    # Lo mismo que la anterior función pero con un solo slider para la velocidad
    phi_0, theta_0, = 3*np.pi/2, 0 # 0, 0
    alpha = np.pi # Parametrizamos con el doble del angulo
    dotphi_0, dottheta_0 = np.sin(alpha*0.5), np.cos(alpha*0.5)

    u = np.linspace(0, 2*np.pi, 100) # phi
    v = np.linspace(0, 2*np.pi, 100) # theta

    U, V = np.meshgrid(u, v)

    X = (b+a*np.sin(U))*np.cos(V)
    Y = (b+a*np.sin(U))*np.sin(V)
    Z = a*np.cos(U)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    superficie = ax.plot_surface(X, Y, Z, alpha=0.5)
    #superficie = ax.plot_wireframe(X, Y, Z)

    # Set an equal aspect ratio
    ax.set_aspect('equal')

    axcolor = 'lightgoldenrodyellow' # definimos el color de los sliders que usaremos más adelante
    ax.margins(x=0)

    sol = resolver(t, a, b, phi_0=phi_0, dotphi_0=dotphi_0, theta_0=theta_0, dottheta_0=dottheta_0)
    phi = sol[:, 0]
    theta = sol[:, 2]

    x = (b+a*np.sin(phi))*np.cos(theta)
    y = (b+a*np.sin(phi))*np.sin(theta)
    z = a*np.cos(phi)

    geodesica = ax.plot(x, y, z, color="red", linewidth=1)

    x_0 = (b+a*np.sin(phi_0))*np.cos(theta_0)
    y_0 = (b+a*np.sin(phi_0))*np.sin(theta_0)
    z_0 = a*np.cos(phi_0)

    punto_inicial = ax.scatter(x_0, y_0, z_0, color="green")

    # ajustamos la posición de la representación para acomodar los sliders
    plt.subplots_adjust(left=0.25, bottom=0.30)

    # hacer un slider vertical para controlar "phi_0". Los cuatro números controlan la posición, el grosor y la altura del slider
    axamp = plt.axes([0.1, 0.25, 0.0225, 0.63], facecolor=axcolor)
    phi_0_slider = Slider(
        ax=axamp,
        label=r"$\phi_0$",
        valmin=0,
        valmax=np.pi*2,
        valinit=phi_0,
        orientation="vertical"
    )

    axamp = plt.axes([0.15, 0.25, 0.0225, 0.63], facecolor=axcolor)
    theta_0_slider = Slider(
        ax=axamp,
        label=r"$\theta_0$",
        valmin=0,
        valmax=np.pi*2,
        valinit=theta_0,
        orientation="vertical"
    )

    axamp = plt.axes([0.20, 0.25, 0.0225, 0.63], facecolor=axcolor)
    alpha_slider = Slider(
        ax=axamp,
        label=r"$\alpha$",
        valmin=0,
        valmax=2*np.pi,
        valinit=alpha,
        orientation="vertical"
    )

    # La función que será llamada cuando se actualice el slider
    def update(val): # para cada vez que se mueve el slider recalculamos la solución
        ax.clear()

        superficie = ax.plot_surface(X, Y, Z, alpha=0.5)

        phi_0, theta_0 = phi_0_slider.val, theta_0_slider.val
        alpha = alpha_slider.val
        dotphi_0, dottheta_0 = np.sin(alpha*0.5), np.cos(alpha*0.5)
        sol = resolver(t, a, b, phi_0=phi_0, dotphi_0=dotphi_0, theta_0=theta_0, dottheta_0=dottheta_0)
        phi = sol[:, 0]
        theta = sol[:, 2]

        x = (b+a*np.sin(phi))*np.cos(theta)
        y = (b+a*np.sin(phi))*np.sin(theta)
        z = a*np.cos(phi)

        geodesica = ax.plot(x, y, z, color="red", linewidth=1)

        x_0 = (b+a*np.sin(phi_0))*np.cos(theta_0)
        y_0 = (b+a*np.sin(phi_0))*np.sin(theta_0)
        z_0 = a*np.cos(phi_0)

        punto_inicial = ax.scatter(x_0, y_0, z_0, color="green")
    
    # Estas lineas registran cuando se cambia el slider
    phi_0_slider.on_changed(update)
    theta_0_slider.on_changed(update)
    alpha_slider.on_changed(update)
    
    # Creamos un boton para resetear los valores del slider
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    def reset(event): # esta función se llamará cuando se pulse el botón
        phi_0_slider.reset()
        theta_0_slider.reset()
        alpha_slider.reset()
    button.on_clicked(reset)

    plt.show()

if __name__ == "__main__":
    # Geodesicas con ángulos constantes
    # Las de theta constante (dottheta = 0)
    # representar_geodesica(t=np.linspace(0, 2*np.pi, 100), a=0.5, b=1, phi_0=0, dotphi_0=1, theta_0=0, dottheta_0=0)

    # dotphi = 0, debe cumplirse que: phi = pi/2, 3*pi/2 para que las geodesicas tengan phi constante
    # Se debe a que cos(phi)(b+asin(phi))dottheta^2 = 0, suponiendo que b > a, solo pasa cuando cos(phi) = 0, phi = (2n+1) * pi/2
    # representar_geodesica(t=np.linspace(0, 10*np.pi, 1000), a=0.5, b=1, phi_0=3*np.pi/2, dotphi_0=0, theta_0=0, dottheta_0=0.5)

    # Otras geodesicas
    # representar_geodesica(t=np.linspace(0, 10*np.pi, 500), a=0.5, b=1, phi_0=np.pi/2, dotphi_0=2**(-1/2), theta_0=0, dottheta_0=2**(-1/2))
    representar_geodesica(t=np.linspace(0, 10*np.pi, 1000), a=0.5, b=1, phi_0=np.pi/2, dotphi_0=1/2, theta_0=0, dottheta_0=(3/4)**0.5)

    # Parametrización de la velocidad como velocidad unitaria
    # alpha = 0.9
    # representar_geodesica(t=np.linspace(0, 10*np.pi, 1000), a=0.5, b=1, phi_0=5*np.pi/4, dotphi_0=np.sin(alpha), theta_0=0, dottheta_0=np.cos(alpha))
    # representar_geodesica_sliders(t=np.linspace(0, 10*np.pi, 1000), a=0.5, b=1)
    # representar_geodesica_sliders2(t=np.linspace(0, 20*np.pi, 2000), a=0.5, b=1)