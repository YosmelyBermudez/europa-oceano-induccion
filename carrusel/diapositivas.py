"""Carrusel de LinkedIn: «¿Hay un océano bajo el hielo de Europa?» (10 diapositivas, 1080×1350)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib.lines import Line2D
from estilo import *
from induccion import europa, esfera_uniforme_analitica, respuesta, respuesta_numerica, R_EUROPA
from campo_inductor import P_SIN

OUT = os.path.join(os.path.dirname(__file__), "png")
os.makedirs(OUT, exist_ok=True)
VJ = dict(lat0=16, lon0=-118)  # vista de Júpiter: parche norte y Gran Mancha Azul
coma = lambda s: s.replace(".", ",")


def guardar(fig, n):
    fig.savefig(f"{OUT}/{n:02d}.png", dpi=DPI, facecolor=BG)
    plt.close(fig)


# ============================================================ 1 · PORTADA
def d01():
    fig, bg = lienzo(1)
    resplandor(bg, 0.40, 0.455, 0.32, alpha=0.55)
    poner_imagen(bg, render_jupiter(900, **VJ), 0.40, 0.455, 0.32)
    # órbita de Europa (arco) y Europa en primer plano
    t = np.linspace(-0.2, np.pi + 0.35, 300)
    bg.plot(0.40 + 0.56 * np.cos(t), 0.455 + 0.12 * np.sin(t) - 0.02, color="#9cc3e0", lw=0.9, alpha=0.55, zorder=2,
            ls=(0, (4, 3)))
    resplandor(bg, 0.82, 0.33, 0.10, color=(0.75, 0.85, 1.0), alpha=0.35)
    poner_imagen(bg, render_europa(520, lat0=5, lon0=40), 0.82, 0.33, 0.10, z=6)
    bg.text(0.82, 0.205, "Europa", color=INK, fontsize=11, ha="center", fontweight="bold", zorder=7)
    bg.text(0.82, 0.183, "(tamaño exagerado)", color=MUTED, fontsize=7.5, ha="center", zorder=7)
    bg.text(0.055, 0.098, "Júpiter: su campo magnético medido por Juno", color=INK2, fontsize=8.5, zorder=7)
    bg.text(0.055, 0.076, "rojo = el campo sale · azul = entra (modelo JRM33)", color=MUTED, fontsize=7.5, zorder=7)
    fig.text(0.055, 0.955, "BITÁCORA DE INVESTIGACIÓN · Nº 1", color=ACCENT, fontsize=9.5, fontweight="bold")
    fig.text(0.055, 0.935, "¿Hay un océano\nbajo el hielo\nde Europa?", color=INK, fontsize=36, fontweight="bold",
             va="top", linespacing=1.05)
    fig.text(0.055, 0.715, "Lo estoy comprobando en Python con la misma física\nque usamos en geofísica para «ver» bajo tierra.",
             color=INK2, fontsize=12.5, va="top", linespacing=1.4)
    fig.text(0.945, 0.07, "Desliza  →", color=ACCENT, fontsize=12, ha="right", fontweight="bold")
    guardar(fig, 1)


# ============================================================ 2 · LA IDEA: MT
def d02():
    fig, bg = lienzo(2)
    cabecera(fig, "La idea", "Europa funciona como un\nsondeo magnetotelúrico gigante",
             "En geofísica medimos cómo responde el subsuelo a un campo\nmagnético que cambia. Aquí hago lo mismo, pero con una luna.")
    ax = fig.add_axes([0.05, 0.30, 0.90, 0.42]); W = 0.90 * 1080 / (0.42 * 1350)
    ax.set_xlim(0, W); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off")
    xm = W / 2
    ax.plot([xm, xm], [0.05, 0.95], color=GRID, lw=0.8)
    ax.text(0.40, 0.98, "EN LA TIERRA", color=ACCENT, fontsize=10, ha="center", fontweight="bold", va="top")
    ax.text(xm + 0.45, 0.98, "EN EUROPA", color=ACCENT, fontsize=10, ha="center", fontweight="bold", va="top")
    def ondas(x0, x1, y):
        xs = np.linspace(x0, x1, 120)
        for k in range(3):
            ax.plot(xs, y - 0.03 * k + 0.012 * np.sin(xs * 38 + k), color="#9cc3e0", lw=1, alpha=0.85)
    # Tierra: fuente
    ax.add_patch(mp.Circle((0.10, 0.82), 0.055, color="#ffcf6b", zorder=3))
    resplandor(ax, 0.10, 0.82, 0.055, color=(1, 0.8, 0.4), alpha=0.5)
    ondas(0.20, 0.78, 0.85)
    ax.text(0.49, 0.72, "el Sol agita la ionosfera y el campo cambia", color=INK2, fontsize=8, ha="center")
    # Tierra: subsuelo
    capas = [(0.58, "#6f7d52"), (0.47, "#86643f"), (0.32, "#2c5f8f"), (0.18, "#4f4036")]
    for y, c in capas:
        ax.add_patch(mp.Rectangle((0.02, 0.04), 0.78, y - 0.04, color=c, zorder=1))
    for wdt in (0.62, 0.42):
        ax.add_patch(mp.Ellipse((0.41, 0.25), wdt, 0.10, fill=False, ec="#ffd27a", lw=1.4, zorder=2, ls="--"))
    ax.text(0.41, 0.25, "roca con agua salada\n(conduce)", color=INK, fontsize=8.3, ha="center", va="center", zorder=3)
    ax.plot([0.41, 0.41], [0.58, 0.65], color=INK, lw=1.5, zorder=3)
    ax.add_patch(mp.Polygon([[0.38, 0.58], [0.44, 0.58], [0.41, 0.63]], color=INK, zorder=3))
    ax.text(0.47, 0.62, "estación\nmagnetotelúrica", color=INK2, fontsize=7.8, va="center")
    ax.text(0.41, 0.10, "corrientes inducidas", color="#ffd27a", fontsize=7.8, ha="center", zorder=3)
    # Europa: fuente
    x0 = xm + 0.06
    resplandor(ax, x0 + 0.06, 0.83, 0.075, alpha=0.4)
    poner_imagen(ax, render_jupiter(300, **VJ), x0 + 0.06, 0.83, 0.075, z=3)
    ondas(x0 + 0.17, W - 0.05, 0.85)
    ax.text(x0 + 0.45, 0.72, "Júpiter gira y el campo cambia", color=INK2, fontsize=8, ha="center")
    # Europa: corte
    cx, cy, r = x0 + 0.38, 0.32, 0.26
    corte_europa(ax, cx, cy, r, ang=(-20, 70), etiquetas=False)
    for wdt in (2 * 0.95 * r, 2 * 0.7 * r):
        ax.add_patch(mp.Ellipse((cx, cy), wdt, 0.07, fill=False, ec="#ffd27a", lw=1.4, zorder=6, ls="--"))
    ax.text(cx, 0.015, "océano salado bajo el hielo (conduce)", color=INK, fontsize=8.3, ha="center")
    # sonda
    sx, sy = W - 0.10, 0.56
    ax.add_patch(mp.Rectangle((sx - 0.02, sy - 0.015), 0.04, 0.03, color="#c9ccd3", zorder=7))
    ax.plot([sx - 0.07, sx + 0.07], [sy, sy], color="#c9ccd3", lw=2.2, zorder=7)
    ax.text(sx, sy + 0.05, "sonda con\nmagnetómetro", color=INK2, fontsize=7.8, ha="center")
    # --- tabla de equivalencias
    filas = [("Fuente del campo", "Sol + ionosfera", "Júpiter"),
             ("Lo que conduce", "Roca con agua salada", "Océano salado"),
             ("Dónde se mide", "Estación en el suelo", "Magnetómetro de la sonda"),
             ("Lo que se deduce", "Conductividad del subsuelo", "Espesor y conductividad del océano")]
    y0 = 0.255
    fig.text(0.29, y0, "TIERRA", color=MUTED, fontsize=8, fontweight="bold")
    fig.text(0.62, y0, "EUROPA", color=MUTED, fontsize=8, fontweight="bold")
    for k, (a, b, c) in enumerate(filas):
        y = y0 - 0.035 - k * 0.038
        fig.text(0.055, y, a, color=INK2, fontsize=9.5)
        fig.text(0.29, y, b, color=INK, fontsize=9.5)
        fig.text(0.62, y, c, color=INK, fontsize=9.5, fontweight="bold")
        fig.add_artist(Line2D([0.055, 0.945], [y - 0.012, y - 0.012], color=GRID, lw=0.6))
    nota(fig, "Esquemas ilustrativos, no a escala.")
    guardar(fig, 2)


# ============================================================ 3 · EL MOTOR
def d03():
    fig, bg = lienzo(3)
    cabecera(fig, "El motor", "Júpiter es un imán inclinado\nque gira cada 9,9 horas",
             "Su ecuador magnético está inclinado ~10°. Europa orbita más lento,\nasí que pasa por encima y por debajo de él cada 11,23 horas.")
    ax = fig.add_axes([0.0, 0.40, 1.0, 0.36]); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-0.75, 0.75); ax.axis("off")
    ax.set_aspect("equal")
    el = np.radians(14)
    def proj(x, y, z):
        return x, z * np.cos(el) + y * np.sin(el), y * np.cos(el) - z * np.sin(el)
    # lámina de plasma (ecuador magnético) inclinada 10° respecto al ecuador de giro
    tilt = np.radians(10)
    t = np.linspace(0, 2 * np.pi, 400)
    def anillo(R):
        x, y = R * np.cos(t), R * np.sin(t)
        return proj(x * np.cos(tilt), y, -x * np.sin(tilt))
    for R, a in [(1.45, 0.10), (1.15, 0.12), (0.85, 0.14), (0.60, 0.16)]:
        X, Y, D = anillo(R)
        ax.fill(X, Y, color="#d98a4a", alpha=a * 0.6, lw=0, zorder=1)
    X, Y, D = anillo(1.45); ax.plot(X, Y, color="#f2a65a", lw=0.8, alpha=0.6, zorder=1)
    # órbita de Europa (ecuatorial)
    Xo, Yo, Do = proj(1.30 * np.cos(t), 1.30 * np.sin(t), 0 * t)
    atras = Do > 0
    ax.plot(np.where(atras, Xo, np.nan), np.where(atras, Yo, np.nan), color="#9cc3e0", lw=1.0, ls=(0, (4, 3)), zorder=2)
    resplandor(ax, 0, 0, 0.42, alpha=0.5)
    poner_imagen(ax, render_jupiter(700, lat0=14, lon0=-118, roll=0), 0, 0, 0.42, z=3)
    ax.plot(np.where(~atras, Xo, np.nan), np.where(~atras, Yo, np.nan), color="#9cc3e0", lw=1.2, ls=(0, (4, 3)), zorder=4)
    # Europa (a la derecha, por encima de la lámina)
    ang = np.radians(-12)
    ex, ey, ed = proj(1.30 * np.cos(ang), 1.30 * np.sin(ang), 0)
    lx, ly, _ = proj(1.30 * np.cos(ang) * np.cos(tilt), 1.30 * np.sin(ang), -1.30 * np.cos(ang) * np.sin(tilt))
    ax.plot([ex, lx], [ey, ly], color=INK2, lw=0.8, ls=":", zorder=5)
    poner_imagen(ax, render_europa(200), ex, ey, 0.07, z=6)
    ax.text(ex, ey + 0.13, "Europa: ahora\npor encima de la lámina", color=INK, fontsize=8.5, ha="center", zorder=7)
    # flecha de giro de la lámina
    tg = np.linspace(np.radians(200), np.radians(250), 50)
    Xg, Yg, _ = proj(1.45 * np.cos(tg) * np.cos(tilt), 1.45 * np.sin(tg), -1.45 * np.cos(tg) * np.sin(tilt))
    ax.annotate("", (Xg[-1], Yg[-1]), (Xg[-8], Yg[-8]), arrowprops=dict(arrowstyle="-|>", color="#f2a65a", lw=1.6))
    ax.plot(Xg, Yg, color="#f2a65a", lw=1.6)
    ax.text(-1.55, 0.50, "La lámina gira con Júpiter.\nEn 5,6 h su lado alto llega\na Europa, que queda debajo.",
            color=INK, fontsize=8.5, zorder=7, va="center")
    ax.text(-1.55, -0.66, "lámina de plasma\n(ecuador magnético, inclinado ~10°)", color="#f2a65a", fontsize=8, zorder=7)
    ax.text(1.55, -0.66, "órbita de Europa\n(no a escala)", color="#9cc3e0", fontsize=8, ha="right", zorder=7)
    # --- curva del campo
    d = np.load(os.path.join(os.path.dirname(__file__), "..", "data", "B_europa_2per.npy"))
    tt, Bx, By, Bz = d
    a2 = fig.add_axes([0.13, 0.15, 0.81, 0.21]); eje_oscuro(a2)
    a2.set_xticks([0, 4, 8, 12, 16, 20])
    a2.plot(tt, By, color=NARANJA, lw=2.4)
    a2.plot(tt, Bx, color=AZUL, lw=2.0)
    a2.axhline(0, color="#3a3f4b", lw=0.8)
    a2.set_xlim(0, tt.max()); a2.set_ylim(-260, 300)
    a2.set_xlabel("Horas"); a2.set_ylabel("Campo que siente Europa (nT)")
    a2.text(0.3, 240, "radial (hacia Júpiter): ±215 nT", color=NARANJA, fontsize=9, fontweight="bold")
    a2.text(0.4, -195, "en la dirección del flujo: ±74 nT", color=AZUL, fontsize=9, fontweight="bold")
    i = np.argmax(By[:120]); a2.annotate("", (tt[i] + P_SIN, 212), (tt[i], 212), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1))
    a2.text(tt[i] + P_SIN / 2, 140, "11,23 h", color=INK, fontsize=9.5, ha="center", fontweight="bold")
    fig.text(0.055, 0.378, "Esa oscilación es la «señal de sondeo»: lo que hace circular corrientes en el océano.",
             color=INK2, fontsize=9.6)
    nota(fig, "Campo calculado con JRM33 (Connerney et al. 2022) + magnetodisco Con2020 (Connerney et al. 2020) en la órbita de Europa.\n"
         "La literatura cita ~200 nT para esta señal (Biersteker et al. 2023).", y=0.048)
    guardar(fig, 3)


# ============================================================ 4 · GALILEO
def d04():
    fig, bg = lienzo(4)
    cabecera(fig, "La pista de Galileo (1996–2000)", "El «imán» de Europa se\ninvierte con el de Júpiter",
             "Un imán propio no cambia de dirección. Uno inducido sí: sigue al campo\nque lo produce. Eso es lo que midió la sonda Galileo.")
    for k, (x0, Bdir, titulo, sub) in enumerate([
            (0.03, +1, "PASO A", "el campo de Júpiter apunta hacia afuera"),
            (0.52, -1, "PASO B", "el campo de Júpiter apunta hacia adentro")]):
        ax = fig.add_axes([x0, 0.29, 0.45, 0.45]); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.7, 1.9); ax.axis("off")
        ax.set_aspect("equal")
        ax.text(0, 1.82, titulo, color=ACCENT, fontsize=10.5, ha="center", fontweight="bold", va="top")
        ax.text(0, 1.55, sub, color=INK2, fontsize=8.6, ha="center", va="top")
        # campo primario: flechas horizontales
        for y in (-1.3, -0.65, 0.65, 1.15):
            xs = (-1.55, 1.55) if Bdir > 0 else (1.55, -1.55)
            ax.annotate("", (xs[1], y), (xs[0], y), arrowprops=dict(arrowstyle="-|>", color="#5f7da3", lw=1.2, alpha=0.8))
        # dipolo inducido opuesto al primario (en la dirección x)
        ang = 180 if Bdir > 0 else 0
        lineas_dipolo(ax, 0, 0, 0.62, ang_m=ang + 90 * 0, Ls=(1.35, 1.9), color="#ffd27a", lw=1.4, z=4)
        resplandor(ax, 0, 0, 0.62, color=(0.75, 0.85, 1.0), alpha=0.25)
        poner_imagen(ax, render_europa(360, lon0=30 + 90 * k), 0, 0, 0.62, z=5)
        ax.annotate("", (0.45 * np.cos(np.radians(ang)), 0), (-0.45 * np.cos(np.radians(ang)), 0),
                    arrowprops=dict(arrowstyle="-|>", color="#ffd27a", lw=3), zorder=6)
        ax.text(0, -1.62, "dipolo inducido →" if ang == 0 else "← dipolo inducido", color="#ffd27a", fontsize=8.8,
                ha="center", fontweight="bold")
    fig.add_artist(Line2D([0.055, 0.15], [0.278, 0.278], color="#5f7da3", lw=1.4))
    fig.text(0.16, 0.272, "campo de Júpiter", color=INK2, fontsize=8.5)
    fig.add_artist(Line2D([0.42, 0.51], [0.278, 0.278], color="#ffd27a", lw=1.4))
    fig.text(0.52, 0.272, "campo inducido por Europa", color=INK2, fontsize=8.5)
    ax = caja(fig, 0.055, 0.075, 0.89, 0.18)
    fig.text(0.08, 0.228, "Por qué importa", color=ACCENT, fontsize=10, fontweight="bold")
    fig.text(0.08, 0.207,
             "Khurana et al. (1998, Nature) interpretaron las perturbaciones cerca de Europa como\n"
             "campos inducidos por una capa conductora bajo la superficie. En enero de 2000,\n"
             "un paso con el campo de Júpiter en la orientación contraria confirmó que el dipolo\n"
             "era inducido y no permanente (Kivelson et al. 2000, Science). La explicación más\n"
             "probable: un océano de agua salada.", color=INK, fontsize=9.0, va="top", linespacing=1.42)
    nota(fig, "Esquema ilustrativo: direcciones simplificadas, no son los datos de los sobrevuelos.", y=0.048)
    guardar(fig, 4)


# ============================================================ 5 · EL MODELO
def d05():
    fig, bg = lienzo(5)
    cabecera(fig, "Mi modelo", "Europa por capas:\nqué responde y cuánto",
             "Programé en Python la respuesta electromagnética de una esfera\ncon hielo, océano y roca. El océano es la única capa que conduce.")
    ax = fig.add_axes([0.02, 0.335, 0.56, 0.40]); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.axis("off")
    ax.set_aspect("equal")
    resplandor(ax, 0, 0, 0.95, color=(0.75, 0.85, 1.0), alpha=0.22)
    pts = corte_europa(ax, 0, 0, 0.95, ang=(-10, 80))
    desc = {"Hielo": ("Corteza de hielo", "0–25 km · no conduce", INK),
            "Océano salado": ("Océano salado", "espesor h · conductividad σ", "#7fb6f0"),
            "Manto rocoso": ("Manto rocoso", "no conduce (en este modelo)", INK),
            "Núcleo de hierro": ("Núcleo metálico", "muy profundo: casi no influye", INK)}
    fig.canvas.draw()
    to_fig = lambda x, y: fig.transFigure.inverted().transform(ax.transData.transform((x, y)))
    for k, (nombre, px, py) in enumerate(pts):
        t1, t2, c = desc[nombre]
        ly = 0.735 - k * 0.088
        fx, fy = to_fig(px, py)
        fig.add_artist(Line2D([fx, 0.60], [fy, ly - 0.004], color=INK2, lw=0.7, alpha=0.8))
        fig.add_artist(mp.Circle((fx, fy), 0.004, transform=fig.transFigure, color=INK))
        fig.text(0.61, ly + 0.004, t1, color=c, fontsize=10.5, fontweight="bold", va="bottom")
        fig.text(0.61, ly - 0.004, t2, color=INK2, fontsize=8.6, va="top")
    fig.text(0.61, 0.375, "Espesores no a escala", color=MUTED, fontsize=7.5)
    # ecuaciones
    ax = caja(fig, 0.055, 0.105, 0.89, 0.215)
    fig.text(0.08, 0.296, "La física (la misma de la magnetotelúrica)", color=ACCENT, fontsize=10, fontweight="bold")
    fig.text(0.08, 0.248, r"$\nabla^2 \mathbf{B} = \mu_0\,\sigma\,\partial \mathbf{B}/\partial t$", color=INK, fontsize=16)
    fig.text(0.50, 0.255, "El campo se difunde en el conductor:\nentra menos cuanto más conduce.", color=INK2,
             fontsize=8.8, va="center", linespacing=1.35)
    fig.text(0.08, 0.180, r"$Q = A\,e^{i\varphi}$", color=INK, fontsize=16)
    fig.text(0.50, 0.185, "A (0 a 1): cuánto campo «expulsa» Europa.\n"
             "φ (0° a 90°): cuánto se retrasa la respuesta.", color=INK2, fontsize=8.8, va="center", linespacing=1.35)
    fig.text(0.08, 0.125, "Se resuelve con funciones esféricas de Bessel capa por capa. Cada par (h, σ) da un A y un φ.",
             color=MUTED, fontsize=8.3)
    nota(fig, "Método clásico de inducción en esferas (p. ej. Zimmer et al. 2000; Biersteker et al. 2023).", y=0.05)
    guardar(fig, 5)


# ============================================================ 6 · VALIDACIÓN
def d06():
    fig, bg = lienzo(6)
    cabecera(fig, "Control de calidad", "Antes de creerle al código,\nlo puse a prueba",
             "Una regla del proyecto: nada se usa sin reproducir primero\nun resultado conocido.")
    sig = np.logspace(-3, 2, 60); sn = sig[::6]
    Qa = np.array([esfera_uniforme_analitica(R_EUROPA, s, P_SIN) for s in sig])
    Qb = np.array([respuesta([R_EUROPA], [s], P_SIN) for s in sig[::3]])
    Qn = np.array([respuesta_numerica([R_EUROPA], [s], P_SIN) for s in sn])
    a = fig.add_axes([0.12, 0.42, 0.82, 0.32]); eje_oscuro(a)
    a.plot(sig, np.abs(Qa), color=AZUL, lw=2.4)
    a.plot(sig[::3], np.abs(Qb), "o", ms=6, mfc="none", mec=NARANJA, mew=1.6)
    a.plot(sn, np.abs(Qn), "x", ms=9, color=INK, mew=1.6)
    a.set_xscale("log"); a.set_xlabel("Conductividad σ (S/m)"); a.set_ylabel("Amplitud A")
    a.set_ylim(-0.02, 1.05)
    a.legend(handles=[Line2D([], [], color=AZUL, lw=2.4, label="Solución analítica clásica"),
                      Line2D([], [], ls="", marker="o", mfc="none", mec=NARANJA, mew=1.6, label="Mi código (capas)"),
                      Line2D([], [], ls="", marker="x", color=INK, mew=1.6, label="Integración numérica")],
             loc="upper left", frameon=False, fontsize=8.6, labelcolor=INK2)
    a.annotate("conductor perfecto:\nexpulsa todo el campo", (60, 0.99), (8, 0.62), color=MUTED, fontsize=8,
               arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
    a.annotate("mal conductor:\ncasi no responde", (1.15e-3, 0.05), (1.1e-3, 0.42), color=MUTED, fontsize=8,
               arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
    # tarjetas de resultados
    tarjetas = [("< 10⁻¹¹", "diferencia entre los\ntres métodos"),
                ("9,925 h", "periodo del marco de los\ndatos de Juno = System III"),
                ("215 nT", "señal en Europa\n(literatura: ~200 nT)")]
    for k, (n, t) in enumerate(tarjetas):
        x = 0.055 + k * 0.302
        caja(fig, x, 0.135, 0.285, 0.20)
        fig.text(x + 0.02, 0.29, "✓", color=AQUA, fontsize=14, fontweight="bold")
        fig.text(x + 0.02, 0.245, n, color=INK, fontsize=19, fontweight="bold")
        fig.text(x + 0.02, 0.215, t, color=INK2, fontsize=8.4, va="top", linespacing=1.35)
    nota(fig, "Datos reales: Juno FGM, 7 jul 2016 (PDS PPI). El modelo JRM33 + magnetodisco reproduce el campo medido\n"
         "(correlación ~0,96 en las componentes horizontales).", y=0.048)
    guardar(fig, 6)


# ============================================================ 7 · RESULTADO
def d07():
    fig, bg = lienzo(7)
    cabecera(fig, "Primer resultado", "Reproduje la cota de\nZimmer et al. (2000)",
             "Galileo exige que Europa responda con A ≥ 0,7. Recorrí miles de océanos\nposibles para ver cuáles cumplen.")
    hs = np.linspace(5, 600, 80); ss = np.logspace(-3, np.log10(30), 80)
    Q = np.array([[europa(h, s, 0, P_SIN) for s in ss] for h in hs]); A = np.abs(Q)
    a = fig.add_axes([0.12, 0.31, 0.70, 0.40]); eje_oscuro(a); a.grid(False)
    cm = LinearSegmentedColormap.from_list("az", ["#0a0f1c", "#132a4d", "#1d4f8f", "#3987e5", "#9cc8f5", "#e8f2fd"])
    cf = a.contourf(ss, hs, A, levels=np.linspace(0, 1, 11), cmap=cm)
    a.contour(ss, hs, A, levels=[0.7], colors=NARANJA, linewidths=3)
    a.axvspan(2.75, 5, color=AQUA, alpha=0.25, lw=0)
    a.text(3.7, 585, "agua de\nmar", color="#0a3b2a", fontsize=8.5, ha="center", va="top", fontweight="bold")
    a.plot(0.0591, 380, "o", ms=11, color=NARANJA, mec="white", mew=1.8, zorder=5)
    a.annotate("mínimo:\nσ = 0,059 S/m", (0.0591, 380), (0.0013, 470), color=INK, fontsize=9.5, fontweight="bold",
               arrowprops=dict(arrowstyle="->", color=INK, lw=1))
    a.text(0.35, 70, "A = 0,7", color=NARANJA, fontsize=11, fontweight="bold")
    a.text(0.0012, 150, "incompatibles\ncon Galileo", color=INK2, fontsize=9)
    a.text(0.6, 300, "compatibles", color="#0a0f1c", fontsize=10, fontweight="bold")
    a.set_xscale("log"); a.set_xlabel("Conductividad del océano σ (S/m)"); a.set_ylabel("Espesor del océano (km)")
    cax = fig.add_axes([0.85, 0.31, 0.025, 0.40]); cb = fig.colorbar(cf, cax=cax)
    cb.set_label("Amplitud A", color=INK2, fontsize=9); cb.ax.tick_params(colors=MUTED, labelsize=8); cb.outline.set_edgecolor(GRID)
    caja(fig, 0.055, 0.105, 0.89, 0.14)
    fig.text(0.08, 0.190, "0,059 S/m", color=NARANJA, fontsize=24, fontweight="bold")
    fig.text(0.08, 0.163, "mi conductividad mínima", color=INK2, fontsize=9)
    fig.text(0.52, 0.190, "0,058 S/m", color=INK, fontsize=24, fontweight="bold")
    fig.text(0.52, 0.163, "Zimmer et al. (2000), Icarus", color=INK2, fontsize=9)
    fig.text(0.08, 0.125, "Diferencia: 2 %. El agua de mar terrestre (2,75–5 S/m) cae muy dentro de la zona compatible.",
             color=INK, fontsize=8.8)
    nota(fig, "Modelo de 3 capas sin corteza de hielo (caso de Zimmer), periodo 11,23 h. La cota 0,058 S/m está citada en Hand y Chyba (2007).", y=0.05)
    guardar(fig, 7)


# ============================================================ 8 · SOBREVUELO
def d08():
    fig, bg = lienzo(8)
    cabecera(fig, "¿Qué vería una sonda?", "Dos océanos distintos dejan\nhuellas distintas",
             "Simulé el campo inducido a lo largo de un paso a 700 km de altura,\ncomo los de Galileo. El tipo de océano cambia la señal en decenas de nT.")
    # esquema de la trayectoria
    ax = fig.add_axes([0.05, 0.505, 0.90, 0.25]); ax.set_xlim(-3.2, 3.2); ax.set_ylim(-1.15, 1.1); ax.axis("off")
    ax.set_aspect("equal")
    lineas_dipolo(ax, 0, 0, 0.75, ang_m=90, Ls=(1.4, 2.0), color="#ffd27a", lw=1.0, alpha=0.6, z=2)
    resplandor(ax, 0, 0, 0.75, color=(0.75, 0.85, 1.0), alpha=0.25)
    poner_imagen(ax, render_europa(380, lon0=70), 0, 0, 0.75, z=4)
    yt = -1.07
    ax.annotate("", (3.0, yt), (-3.0, yt), arrowprops=dict(arrowstyle="-|>", color="#c9ccd3", lw=1.6), zorder=5)
    ax.text(-3.0, yt + 0.12, "trayectoria simulada", color=INK2, fontsize=8.5)
    ax.plot(0, yt, "o", color=ACCENT, ms=6, zorder=6); ax.text(0.12, yt + 0.12, "máxima aproximación (700 km)", color=ACCENT, fontsize=8)
    ax.annotate("", (2.6, 0.95), (2.6, 0.35), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.2))
    ax.text(2.45, 0.65, "hacia\nJúpiter", color=INK2, fontsize=8, ha="right", va="center")
    ax.text(-1.55, 0.75, "campo inducido\n(dipolo)", color="#ffd27a", fontsize=8, ha="right", va="center")
    # curvas
    d = np.load(os.path.join(os.path.dirname(__file__), "..", "data", "B_europa_2per.npy"))
    tt, Bx, By, Bz = d; N = len(tt)
    F = np.fft.rfft(np.vstack([Bx, By, Bz]), axis=1); Bhat = F[:, 2] * 2 / N
    w = 2 * np.pi / P_SIN; t0 = tt[np.argmax(By)]
    R = R_EUROPA / 1e3; tm = np.linspace(-15, 15, 400); s = tm * 60 * 6.0
    pos = np.stack([s, -(R + 700) * np.ones_like(s), 0 * s]); rn = np.linalg.norm(pos, axis=0); rh = pos / rn
    a = fig.add_axes([0.12, 0.17, 0.82, 0.30]); eje_oscuro(a)
    casos = [("Océano mínimo · σ = 0,06 S/m, 380 km", europa(380, 0.06, 0, P_SIN), AZUL),
             ("Agua de mar · σ = 3 S/m, 100 km, bajo 25 km de hielo", europa(100, 3.0, 25, P_SIN), NARANJA)]
    for lab, Q, c in casos:
        M = -0.5 * np.real(Q * Bhat * np.exp(1j * w * t0))
        Bi = (R / rn) ** 3 * (3 * (M[:, None] * rh).sum(0) * rh - M[:, None])
        a.plot(tm, Bi[1], color=c, lw=2.6)
    a.axhline(0, color="#3a3f4b", lw=0.8)
    a.text(-14.5, -64, "Sin océano: línea plana en 0", color=MUTED, fontsize=8.5)
    a.text(1.5, -60, "agua de mar\n(A = 0,90)", color=NARANJA, fontsize=9.5, fontweight="bold")
    a.text(3.0, -33, "océano mínimo (A = 0,70)", color=AZUL, fontsize=9.5, fontweight="bold")
    a.annotate("", (0, -62.5), (0, -39), arrowprops=dict(arrowstyle="<->", color=INK, lw=1))
    a.text(-0.5, -51, "~24 nT", color=INK, fontsize=10, fontweight="bold", ha="right")
    a.set_xlabel("Minutos desde la máxima aproximación"); a.set_ylabel("Campo inducido, componente\nhacia Júpiter (nT)")
    a.set_ylim(-70, 10)
    nota(fig, "Simulación ilustrativa (no es un sobrevuelo real). Solo el dipolo inducido: en los datos reales también hay\n"
         "perturbaciones del plasma del mismo orden, que habrá que modelar (Schilling et al. 2007).", y=0.048)
    guardar(fig, 8)


# ============================================================ 9 · REFERENCIAS
def d09():
    fig, bg = lienzo(9)
    cabecera(fig, "En qué me apoyo", "Papers, datos y herramientas",
             "Todo el trabajo parte de resultados publicados y datos abiertos de la NASA.")
    grupos = [
        ("EL OCÉANO DE EUROPA", [
            ("Khurana et al. (1998)", "Nature 395: campos inducidos como evidencia de océanos en Europa y Calisto."),
            ("Kivelson et al. (2000)", "Science: el dipolo se invierte con el campo de Júpiter."),
            ("Zimmer, Khurana y Kivelson (2000)", "Icarus 147: cotas de conductividad del océano."),
            ("Schilling, Neubauer y Saur (2007)", "Icarus: efecto del plasma; σ × espesor ≥ 50 S/m·km."),
            ("Hand y Chyba (2007)", "Icarus 189: salinidad del océano a partir de la inducción.")]),
        ("MÉTODO Y FUTURO", [
            ("Biersteker et al. (2023)", "Planetary Science Journal: inversión bayesiana para Europa Clipper.")]),
        ("CAMPO DE JÚPITER", [
            ("Connerney et al. (2022)", "JGR Planets: modelo JRM33 con 33 órbitas de Juno."),
            ("Connerney et al. (2020)", "JGR Space Physics: modelo del magnetodisco Con2020.")]),
        ("DATOS Y CÓDIGO", [
            ("NASA PDS PPI", "Magnetómetro de Juno (JNO-J-3-FGM-CAL-V1.0) y, próximamente, Galileo MAG."),
            ("Python", "NumPy, SciPy, matplotlib y JupiterMag (modelos JRM33 y Con2020).")]),
    ]
    y = 0.80
    for titulo, items in grupos:
        fig.text(0.055, y, titulo, color=ACCENT, fontsize=9.2, fontweight="bold"); y -= 0.034
        for a, b in items:
            fig.text(0.075, y, a, color=INK, fontsize=10, fontweight="bold")
            fig.text(0.075, y - 0.022, b, color=INK2, fontsize=8.8)
            y -= 0.056
        y -= 0.012
    nota(fig, "Referencias completas en el README del proyecto.", y=0.05)
    guardar(fig, 9)


# ============================================================ 10 · LO QUE VIENE
def d10():
    fig, bg = lienzo(10)
    resplandor(bg, 0.80, 0.30, 0.20, alpha=0.4)
    poner_imagen(bg, render_jupiter(500, **VJ), 0.80, 0.30, 0.20)
    poner_imagen(bg, render_europa(260, lon0=200), 0.52, 0.18, 0.05, z=6)
    cabecera(fig, "Lo que viene", "Del modelo a los datos reales",
             "La entrada 1 está lista: modelo construido y validado.\nEstos son los siguientes pasos:")
    pasos = [("1", "Datos reales de Galileo", "Sobrevuelos E4, E14 y E26 del archivo PDS."),
             ("2", "Ajustar A y φ", "Restar el campo de Júpiter y ajustar la respuesta a las mediciones."),
             ("3", "Repetir la prueba de 2000", "Comparar pasos con el campo de Júpiter en sentidos opuestos."),
             ("4", "Mirar al futuro", "Europa Clipper (NASA) y JUICE (ESA) medirán esto con más precisión.")]
    for k, (n, t, d) in enumerate(pasos):
        y = 0.75 - k * 0.085
        fig.add_artist(mp.Ellipse((0.085, y + 0.005), 0.044, 0.044 * 1080 / 1350, transform=fig.transFigure, color=ACCENT, zorder=3))
        fig.text(0.085, y + 0.005, n, color=BG, fontsize=11, fontweight="bold", ha="center", va="center", zorder=4)
        fig.text(0.13, y + 0.013, t, color=INK, fontsize=11.5, fontweight="bold", va="center")
        fig.text(0.13, y - 0.014, d, color=INK2, fontsize=9, va="center")
    fig.text(0.055, 0.34, "Iré contando cada paso\nen esta bitácora.", color=INK, fontsize=17, fontweight="bold",
             va="top", linespacing=1.15)
    fig.text(0.055, 0.25, "¿Trabajas en ciencias planetarias,\ngeofísica o magnetismo espacial?\nMe encantaría conversar.",
             color=INK2, fontsize=11, va="top", linespacing=1.4)
    guardar(fig, 10)


if __name__ == "__main__":
    which = [int(a) for a in sys.argv[1:]] or list(range(1, 11))
    for n in which:
        globals()[f"d{n:02d}"]()
        print("ok", n)
