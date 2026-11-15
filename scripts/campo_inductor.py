import numpy as np, JupiterMag as jm
A_EUR = 9.38          # RJ (radio orbital de Europa ≈ 671 100 km / 71 492 km)
P_JUP = 9.9250        # h, System III
P_EUR = 85.228        # h, periodo orbital
P_SIN = 1/(1/P_JUP - 1/P_EUR)

def campo_en_europa(t_h, lon0=0.0, disco=True):
    """B (nT) que siente Europa en coordenadas EPhiO: x = flujo de corrotación,
    y = hacia Júpiter, z = eje de giro de Júpiter."""
    lam = np.radians(lon0) - 2*np.pi*t_h/P_SIN     # longitud este en System III
    x, y, z = A_EUR*np.cos(lam), A_EUR*np.sin(lam), 0*lam
    jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=True, Degree=18)
    B = np.array(jm.Internal.Field(x, y, z))
    if disco:
        B = B + np.array(jm.Con2020.Field(x, y, z))
    rhat = np.array([np.cos(lam), np.sin(lam), 0*lam]); phat = np.array([-np.sin(lam), np.cos(lam), 0*lam])
    return np.array([(B*phat).sum(0), -(B*rhat).sum(0), B[2]]), lam
