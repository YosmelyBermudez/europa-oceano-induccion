"""Respuesta de inducción electromagnética (grado n=1) de una esfera de capas.

Convención: dependencia temporal e^{iωt}. Q = A e^{-iφ} es la razón entre el
dipolo inducido y el campo primario uniforme, normalizada para que un
conductor perfecto en la superficie dé A = 1, φ = 0 (convención de Zimmer et al. 2000).
El campo inducido fuera de la luna es
    B_ind = -(A e^{-iφ}/2) (R/r)^3 [3 (B_p·r̂) r̂ - B_p].

Método: en cada capa conductora el escalar poloidal S(r) es combinación de las
funciones esféricas de Bessel modificadas i_1(kr), k_1(kr), con k = sqrt(iωμ0σ);
en capas aislantes, de r y r^-2. Se imponen S y dS/dr continuos en cada interfaz.
"""
import numpy as np
from scipy.special import spherical_in, spherical_kn
from scipy.integrate import solve_ivp

MU0 = 4e-7 * np.pi
R_EUROPA = 1560.8e3  # m


def _base(r, sigma, omega):
    """Dos soluciones independientes y sus derivadas en r para una capa."""
    if sigma == 0:
        return np.array([r, r**-2.0]), np.array([1.0, -2.0 * r**-3.0])
    k = np.sqrt(1j * omega * MU0 * sigma)
    z = k * r
    f = np.array([spherical_in(1, z), spherical_kn(1, z)])
    df = k * np.array([spherical_in(1, z, derivative=True),
                       spherical_kn(1, z, derivative=True)])
    return f, df


def respuesta(radios, sigmas, periodo_h):
    """radios: radios externos de cada capa (m), de dentro hacia fuera; el último
    es la superficie. sigmas: conductividad de cada capa (S/m).
    Devuelve Q = A e^{-iφ} (complejo)."""
    omega = 2 * np.pi / (periodo_h * 3600.0)
    radios = np.asarray(radios, float)
    # capa central: solo la solución regular
    c = np.array([1.0, 0.0], complex)
    for j in range(len(radios) - 1):
        rb = radios[j]
        f, df = _base(rb, sigmas[j], omega)
        S, dS = c @ f, c @ df
        f2, df2 = _base(rb, sigmas[j + 1], omega)
        M = np.array([f2, df2])
        c = np.linalg.solve(M, np.array([S, dS]))
        # reescalado para evitar desbordes (la respuesta es una razón)
        c = c / np.max(np.abs(c))
    R = radios[-1]
    f, df = _base(R, sigmas[-1], omega)
    S, dS = c @ f, c @ df
    # fuera: S = a r + g r^-2
    a, g = np.linalg.solve(np.array([[R, R**-2], [1.0, -2 * R**-3]]), [S, dS])
    return -g / (a * R**3)


def europa(h_oceano_km, sigma_oceano, hielo_km=0.0, periodo_h=11.23,
           sigma_manto=0.0):
    """Modelo de 3 capas: interior (manto/núcleo) | océano | hielo aislante."""
    R = R_EUROPA
    r_top = R - hielo_km * 1e3
    r_bot = r_top - h_oceano_km * 1e3
    radios = [r_bot, r_top, R]
    sig = [sigma_manto, sigma_oceano, 0.0]
    return respuesta(radios, sig, periodo_h)


def A_phi(Q):
    return np.abs(Q), np.degrees(np.angle(Q))


# ---------------------------------------------------------------------------
# Validaciones independientes
# ---------------------------------------------------------------------------
def esfera_uniforme_analitica(a, sigma, periodo_h):
    """Solución cerrada clásica para una esfera uniforme (radio a)."""
    omega = 2 * np.pi / (periodo_h * 3600.0)
    x = np.sqrt(1j * omega * MU0 * sigma) * a
    return 1 + 3 / x**2 - 3 / (x * np.tanh(x))


def respuesta_numerica(radios, sigmas, periodo_h, n_sub=4000):
    """Integra la EDO S'' + 2S'/r - 2S/r^2 = k^2 S capa por capa con solve_ivp."""
    omega = 2 * np.pi / (periodo_h * 3600.0)
    radios = np.asarray(radios, float)
    r0 = radios[0] * 1e-3 if sigmas[0] > 0 else radios[0]
    # arranque: en r pequeño (o en el borde del interior aislante) S ~ r
    y = np.array([r0, 1.0], complex)
    if sigmas[0] > 0:
        k = np.sqrt(1j * omega * MU0 * sigmas[0])
        y = np.array([spherical_in(1, k * r0), k * spherical_in(1, k * r0, True)])
        tramos = [(r0, radios[0], sigmas[0])]
    else:
        tramos = []
    tramos += [(radios[j], radios[j + 1], sigmas[j + 1]) for j in range(len(radios) - 1)]
    for (ra, rb, s) in tramos:
        k2 = 1j * omega * MU0 * s

        def f(r, Y):
            return [Y[1], k2 * Y[0] - 2 * Y[1] / r + 2 * Y[0] / r**2]
        sol = solve_ivp(f, (ra, rb), y, method="DOP853", rtol=1e-10, atol=1e-30)
        y = sol.y[:, -1]
        y = y / np.max(np.abs(y))
    R = radios[-1]
    a, g = np.linalg.solve(np.array([[R, R**-2], [1.0, -2 * R**-3]]), y)
    return -g / (a * R**3)
