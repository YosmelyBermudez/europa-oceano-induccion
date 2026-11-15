"""Lectura de datos de Galileo MAG (EPhiO) y ajuste del dipolo inducido de Europa."""
import numpy as np, pandas as pd, os
DIR = os.path.join(os.path.dirname(__file__), "..", "data", "galileo")
PASOS = {"E4": "ORB04", "E14": "ORB14", "E26": "ORB26"}

def leer(paso):
    f = os.path.join(DIR, f"{PASOS[paso]}_EUR_EPHIO.TAB")
    d = pd.read_csv(f, sep=r"\s+", header=None, names=["t", "Bx", "By", "Bz", "B", "x", "y", "z"])
    d["t"] = pd.to_datetime(d["t"])
    d["r"] = np.sqrt(d.x**2 + d.y**2 + d.z**2)
    ica = d.r.idxmin()
    d["min"] = (d.t - d.t[ica]).dt.total_seconds() / 60
    return d

def fondo(d, rcorte=3.5, grado=3):
    """Campo de fondo: polinomio en el tiempo ajustado solo lejos de Europa (r > rcorte)."""
    m = d.r > rcorte
    out = {}
    for c in ["Bx", "By", "Bz"]:
        p = np.polyfit(d["min"][m], d[c][m], grado)
        out[c] = np.polyval(p, d["min"])
    return pd.DataFrame(out)

def campo_dipolo(M, x, y, z):
    """Campo (nT) de un dipolo centrado; M en nT (campo ecuatorial en la superficie = M)."""
    P = np.vstack([x, y, z]); r = np.linalg.norm(P, axis=0); rh = P / r
    return (1 / r) ** 3 * (3 * (M[:, None] * rh).sum(0) * rh - M[:, None])

def matriz_dipolo(x, y, z):
    cols = [campo_dipolo(e, x, y, z) for e in np.eye(3)]
    return np.stack(cols, axis=-1)   # (3, N, 3)

def ajustar(d, rcorte=3.5, rmax=4.0, grado=3):
    bg = fondo(d, rcorte, grado)
    dB = d[["Bx", "By", "Bz"]].values - bg.values
    m = (d.r < rmax).values
    G = matriz_dipolo(d.x.values[m], d.y.values[m], d.z.values[m])
    Gf = G.transpose(1, 0, 2).reshape(-1, 3); y = dB[m].ravel()
    M, *_ = np.linalg.lstsq(Gf, y, rcond=None)
    ica = d.r.idxmin()
    Bp = np.array([bg.Bx[ica], bg.By[ica], 0.0])         # parte variable del campo primario
    g = (Gf @ (-0.5 * Bp))                                # respuesta con A = 1
    A = float(g @ y / (g @ g))
    res_libre = y - Gf @ M; res_A = y - A * g
    return dict(M=M, A=A, Bp=Bp, bg=bg, dB=dB, mask=m,
                rms_libre=np.sqrt(np.mean(res_libre**2)), rms_A=np.sqrt(np.mean(res_A**2)),
                rms_dato=np.sqrt(np.mean(y**2)), ica=ica)
