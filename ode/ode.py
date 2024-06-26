



def euler(x0,t,f):

    """Metodo Euler para Ecuaciones Diferenciales Ordinarias con condicion inicial

    Example:
        >>>euler(0,t,func)
        [ 0.          0.          1.18623077 -0.15217513 -0.86222182]

    Arguments:

        x0 -- condicion inicial 
        t (array)  -- arreglo de intervalo a evaluar 
        f -- funcion en estudio

    Returns:
        array :  Un arreglo que contiene la solucion de la funcion en estudio evaluada en los puntos dados.
    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size - 1):
        x[i+1] = x[i] + (h * f(x[i],t[i]))
    return x


def rk2(x0,t,f):
    """Metodo Runge-Kutta de 2do orden para Ecuaciones Diferenciales Ordinarias con condicion inicial

    Example:
        >>>rk2(0,t,func)
        [0.         0.         0.94312516 1.20968194 2.17349484]

    Args:
        x0 -- Condicion inicial
        t (array) -- Arreglo del intervalo a evaluar
        func -- Funcion en estudio



    Returns:
        array :  Un arreglo que contiene la solucion de la funcion en estudio evaluada en los puntos dados.

    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i] + k1/2 , t[i] + k1/2)
        x[i+1] = x[i] + k2
    return x


def rk4(x0,t,f):
    """Metodo Runge-Kutta de 4to orden para Ecuaciones Diferenciales Ordinarias con condicion inicial

    Example:
        >>>rk4(0,t,func)
        [ 0.          0.19770513  0.63206696  0.2423863  -0.35739581]

    Args:
        x0 -- Condicion inicial
        t (array) -- Arreglo del intervalo a evaluar
        func -- Funcion en estudio

    Returns:
        array :  Un arreglo que contiene la solucion de la funcion en estudio evaluada en los puntos dados.

    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i] + k1/2 , t[i] + k1/2)
        k3 = h * f(x[i] + k2/2 , t[i] + k2/2)
        k4 = h * f(x[i] + k3 , t[i] + h)
        x[i+1] = x[i] + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x

