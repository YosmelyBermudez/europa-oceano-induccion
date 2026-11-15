"""Estilo «revista de astronomía» para el carrusel: fondo estelar, renders de
Júpiter (Br de JRM33) y de Europa, cortes de capas y utilidades de maquetación."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.patches as mp

for f in font_manager.findSystemFonts():
    if "/inter/" in f.lower():
        font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = "Inter"

# --- paleta -------------------------------------------------------------------
BG = "#05060a"
INK = "#f3f1ec"; INK2 = "#b9bcc6"; MUTED = "#80838d"
ACCENT = "#f2a65a"          # ámbar (titulares de sección)
AZUL, NARANJA, AQUA = "#3987e5", "#d95926", "#199e70"
PANEL = "#0c0f16"; GRID = "#252a35"
W_IN, H_IN, DPI = 7.2, 9.0, 150   # 1080 × 1350 px
TOTAL = 10
AUTORA = "Yosmely Bermúdez · Ingeniera Geofísica · Ciencia de datos"


def lienzo(n, seed=None):
    fig = plt.figure(figsize=(W_IN, H_IN), dpi=DPI, facecolor=BG)
    bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, 1); bg.set_ylim(0, 1.25); bg.axis("off")
    bg.set_facecolor(BG)
    rng = np.random.default_rng(seed if seed is not None else 100 + n)
    k = 420
    x, y = rng.random(k), rng.random(k) * 1.25
    s = rng.exponential(0.35, k).clip(0.02, 3.0)
    col = rng.choice(["#ffffff", "#cfd9ff", "#ffe2c2"], k, p=[0.6, 0.25, 0.15])
    bg.scatter(x, y, s=s, c=col, alpha=0.75, lw=0, zorder=0)
    # pie
    fig.text(0.055, 0.028, AUTORA, color=MUTED, fontsize=7.6)
    fig.text(0.945, 0.028, f"{n:02d} / {TOTAL:02d}", color=MUTED, fontsize=7.6, ha="right")
    return fig, bg


def cabecera(fig, kicker, titulo, sub=None, y=0.955, size=23):
    fig.text(0.055, y, kicker.upper(), color=ACCENT, fontsize=9, fontweight="bold")
    t = fig.text(0.055, y - 0.018, titulo, color=INK, fontsize=size, fontweight="bold",
                 va="top", linespacing=1.12)
    if sub:
        nlin = titulo.count("\n") + 1
        fig.text(0.055, y - 0.03 - nlin * size * 0.0021 - 0.012, sub, color=INK2, fontsize=11.2,
                 va="top", linespacing=1.4)
    return t


def nota(fig, txt, y=0.058):
    fig.text(0.055, y, txt, color=MUTED, fontsize=7.4, va="bottom", linespacing=1.35)


def resplandor(ax, cx, cy, r, color=(1.0, 0.55, 0.2), alpha=0.5, n=40):
    for i in range(n):
        f = 1 + 1.4 * i / n
        ax.add_patch(mp.Circle((cx, cy), r * f, color=color, alpha=alpha * (1 - i / n) ** 2.2 / n * 3, lw=0, clip_on=False,
                               zorder=1))


# --- render ortográfico genérico -------------------------------------------------
def _ortho(N, lat0, lon0, roll=0.0):
    yy, xx = np.mgrid[1:-1:N * 1j, -1:1:N * 1j]
    rr = xx**2 + yy**2
    dentro = rr <= 1
    zz = np.sqrt(np.clip(1 - rr, 0, 1))
    c, s = np.cos(np.radians(roll)), np.sin(np.radians(roll))
    xr, yr = c * xx - s * yy, s * xx + c * yy
    la, lo = np.radians(lat0), np.radians(lon0)
    ev = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    ee = np.array([-np.sin(lo), np.cos(lo), 0.0])
    en = np.cross(ev, ee)
    P = xr[..., None] * ee + yr[..., None] * en + zz[..., None] * ev
    lat = np.degrees(np.arcsin(np.clip(P[..., 2], -1, 1)))
    lon = np.degrees(np.arctan2(P[..., 1], P[..., 0]))
    return xx, yy, zz, dentro, lat, lon, P


def _sombra(xx, yy, zz, luz=(-0.55, 0.45, 0.70), amb=0.16):
    L = np.array(luz) / np.linalg.norm(luz)
    d = np.clip(xx * L[0] + yy * L[1] + zz * L[2], 0, 1)
    return amb + (1 - amb) * d**0.9 * (0.55 + 0.45 * zz**0.35)


_BR = None


def mapa_br(nlat=181, nlon=361):
    global _BR
    if _BR is None:
        import JupiterMag as jm
        jm.Internal.Config(Model="jrm33", CartesianIn=False, CartesianOut=False, Degree=18)
        lat = np.linspace(-90, 90, nlat); lon = np.linspace(-180, 180, nlon)
        LA, LO = np.meshgrid(lat, lon, indexing="ij")
        th = np.radians(90 - LA).ravel(); ph = np.radians(LO).ravel()
        Br = np.array(jm.Internal.Field(np.ones_like(th), th, ph))[0].reshape(LA.shape)
        _BR = (lat, lon, Br / 1e3)  # µT
    return _BR


CMAP_BR = LinearSegmentedColormap.from_list("br", ["#0b2f63", "#2f6fb0", "#9cc3e0", "#efe9e2",
                                                   "#f2b08f", "#c8432f", "#6b0a1f"])


def render_jupiter(N=700, lat0=18, lon0=200, roll=0, vmax=1600, grat=True):
    lat, lon, Br = mapa_br()
    xx, yy, zz, dentro, la, lo, P = _ortho(N, lat0, lon0, roll)
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((lat, lon), Br)
    val = f(np.stack([np.clip(la, -90, 90), np.clip(lo, -180, 180)], -1))
    rgb = CMAP_BR(Normalize(-vmax, vmax)(val))[..., :3]
    # bandas sutiles de nubes para que «parezca» Júpiter
    bandas = 1 + 0.035 * np.sin(np.radians(la) * 14) + 0.02 * np.sin(np.radians(la) * 31 + 1)
    sh = _sombra(xx, yy, zz)[..., None] * bandas[..., None]
    img = np.clip(rgb * sh, 0, 1)
    if grat:
        px = 180.0 / N  # ~grados por píxel en el centro
        g = np.minimum(np.abs(((la + 15) % 30) - 15), np.abs(((lo + 15) % 30) - 15) * np.cos(np.radians(la)))
        wgt = np.clip(1 - g / (px * 0.9), 0, 1) * 0.16
        img = img * (1 - wgt[..., None]) + 0.95 * wgt[..., None]
    a = np.where(dentro, 1.0, 0.0)
    # borde suave
    rr = np.sqrt(xx**2 + yy**2); a = np.clip((1 - rr) * N / 3, 0, 1) * a
    return np.dstack([img, a])


def render_europa(N=500, lat0=10, lon0=30, seed=7):
    xx, yy, zz, dentro, la, lo, P = _ortho(N, lat0, lon0)
    rng = np.random.default_rng(seed)
    base = np.array([0.90, 0.86, 0.79])
    # moteado (manchas pardas)
    mot = np.zeros(la.shape)
    for _ in range(14):
        k = rng.normal(size=3); k /= np.linalg.norm(k)
        mot += np.sin((P @ k) * rng.uniform(3, 9) + rng.uniform(0, 6))
    mot = (mot - mot.min()) / (np.ptp(mot) + 1e-9)
    col = base[None, None] * (1 - 0.18 * mot[..., None]) + np.array([0.55, 0.38, 0.25]) * 0.18 * mot[..., None]
    # lineae: arcos de círculos máximos ondulados
    lin = np.zeros(la.shape)
    for _ in range(26):
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        q = rng.normal(size=3); q /= np.linalg.norm(q)
        d = np.abs(P @ n + 0.03 * np.sin((P @ q) * 20))
        w = rng.uniform(0.004, 0.012)
        lin = np.maximum(lin, np.exp(-(d / w) ** 2) * rng.uniform(0.35, 0.8))
    col = col * (1 - lin[..., None]) + np.array([0.56, 0.30, 0.18]) * lin[..., None]
    img = np.clip(col * _sombra(xx, yy, zz, amb=0.12)[..., None], 0, 1)
    rr = np.sqrt(xx**2 + yy**2); a = np.clip((1 - rr) * N / 3, 0, 1) * dentro
    return np.dstack([img, a])


def poner_imagen(ax, img, cx, cy, r, z=3):
    ax.imshow(img, extent=[cx - r, cx + r, cy - r, cy + r], zorder=z, interpolation="bilinear")


def corte_europa(ax, cx, cy, r, ang=(0, 90), capas=None, etiquetas=True, fs=9.5, texto_lado="derecha"):
    """Europa con una cuña recortada que muestra sus capas (espesores no a escala)."""
    poner_imagen(ax, render_europa(420), cx, cy, r, z=3)
    if capas is None:
        capas = [("Hielo", 1.00, "#dfe8ef"), ("Océano salado", 0.90, "#1f5fa8"),
                 ("Manto rocoso", 0.76, "#8a5a3c"), ("Núcleo de hierro", 0.42, "#9a9ea6")]
    a0, a1 = ang
    ax.add_patch(mp.Wedge((cx, cy), r * 1.005, a0, a1, color=BG, zorder=4, lw=0))
    for nombre, f, c in capas:
        ax.add_patch(mp.Wedge((cx, cy), r * f, a0, a1, color=c, zorder=5, lw=0))
    # brillo del océano
    ax.add_patch(mp.Wedge((cx, cy), r * 0.90, a0, a1, width=r * 0.03, color="#4f95e0", zorder=5.5, lw=0, alpha=0.8))
    puntos = []
    if etiquetas:
        am = np.radians((a0 + a1) / 2)
        prev = [1.0] + [c[1] for c in capas[1:]] + [0.0]
        for (nombre, f, c), fin in zip(capas, prev[1:]):
            rm = (f + fin) / 2 if nombre != "Núcleo de hierro" else f * 0.5
            px, py = cx + r * rm * np.cos(am), cy + r * rm * np.sin(am)
            puntos.append((nombre, px, py))
    return puntos


def lineas_dipolo(ax, cx, cy, r, ang_m, Ls=(1.5, 2.2, 3.2), color="#ffd27a", lw=1.2, alpha=0.9, z=2, rmin=1.02):
    """Líneas de un dipolo con momento en dirección ang_m (grados) centrado en (cx, cy)."""
    a = np.radians(ang_m)
    for L in Ls:
        th = np.linspace(0.02, np.pi - 0.02, 400)
        rr = L * np.sin(th) ** 2
        for lado in (1, -1):
            x = rr * np.cos(th); y = lado * rr * np.sin(th)
            m = rr >= rmin
            X = cx + r * (x * np.cos(a) - y * np.sin(a)); Y = cy + r * (x * np.sin(a) + y * np.cos(a))
            Xm, Ym = np.where(m, X, np.nan), np.where(m, Y, np.nan)
            ax.plot(Xm, Ym, color=color, lw=lw, alpha=alpha, zorder=z, solid_capstyle="round")


def eje_oscuro(ax):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.grid(color=GRID, lw=0.6)
    ax.xaxis.label.set_color(INK2); ax.yaxis.label.set_color(INK2)
    ax.xaxis.label.set_size(9.5); ax.yaxis.label.set_size(9.5)


def caja(fig, x, y, w, h, color=PANEL, borde=GRID, z=0.5):
    ax = fig.add_axes([x, y, w, h]); ax.axis("off")
    ax.add_patch(mp.FancyBboxPatch((0.005, 0.01), 0.99, 0.98, boxstyle="round,pad=0,rounding_size=0.04",
                                   transform=ax.transAxes, fc=color, ec=borde, lw=0.8, alpha=0.92))
    ax.set_zorder(z)
    return ax
