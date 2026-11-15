"""Figuras de la entrada 2: datos reales de Galileo (E4, E14, E26)."""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from galileo import leer, ajustar, campo_dipolo, PASOS

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "lines.linewidth": 1.8})
OUT = "../figs/"
FECHAS = {"E4": "19 dic 1996", "E14": "29 mar 1998", "E26": "3 ene 2000"}


def titulo(fig, t, sub):
    fig.text(0.012, 0.985, t, fontsize=13.5, fontweight="bold", color=INK, va="top")
    fig.text(0.012, 0.935, sub, fontsize=10, color=INK2, va="top")


def nota(fig, txt):
    fig.text(0.012, 0.012, txt, fontsize=8, color=MUTED, va="bottom")


datos = {p: leer(p) for p in PASOS}
fits = {p: ajustar(datos[p], 3.5, 4.0, 3) for p in PASOS}

# robustez: 27 combinaciones de procesamiento
rob = {}
for p in PASOS:
    As = [ajustar(datos[p], rc, rm, g)["A"] for rc, rm, g in itertools.product([3.0, 3.5, 4.0], [2.5, 3.0, 4.0], [2, 3, 4])]
    rob[p] = (np.median(As), min(As), max(As))

# --------------------------------------------------- Fig 6: datos crudos + fondo
fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), sharey=False)
fig.subplots_adjust(top=0.72, bottom=0.18, left=0.06, right=0.99, wspace=0.22)
titulo(fig, "Paso 1: separar el campo de Júpiter de la huella de Europa",
       "Componente By (hacia Júpiter) medida por Galileo. Gris: tramo cerca de Europa (r < 3,5 radios). La curva naranja es el fondo de Júpiter,\n"
       "ajustado solo con los tramos lejanos. Lo que queda entre los datos y esa curva es la perturbación que produce Europa.")
for ax, p in zip(axs, PASOS):
    d, f = datos[p], fits[p]
    cerca = d.r < 3.5
    ax.axvspan(d["min"][cerca].min(), d["min"][cerca].max(), color="#ecebe6", lw=0)
    ax.plot(d["min"], d.By, color=INK, lw=1.0)
    ax.plot(d["min"], f["bg"].By, color=C2, lw=2.2)
    ax.axvline(0, color=MUTED, lw=0.8, ls=":")
    ax.set_title(f"{p} · {FECHAS[p]} · {(d.r.min() - 1) * 1560:.0f} km\nfondo By en la máxima aproximación: {f['Bp'][1]:+.0f} nT".replace("-", "−"), fontsize=10)
    ax.set_xlabel("Minutos desde la máxima aproximación")
axs[0].set_ylabel("By (nT)")
axs[0].legend(handles=[Line2D([], [], color=INK, lw=1, label="Galileo (medido)"),
                       Line2D([], [], color=C2, lw=2.2, label="fondo de Júpiter (ajustado)")],
              loc="upper right", frameon=False, fontsize=8.5)
nota(fig, "Datos: Galileo MAG alta resolución (GO-J-MAG-3-RDR-HIGHRES-V1.0, PDS PPI), 1/3 s, coordenadas EPhiO. Fondo: polinomio cúbico en el tiempo ajustado con r > 3,5 RE.")
fig.savefig(OUT + "fig6_galileo_fondo.png", dpi=160); plt.close(fig)

# --------------------------------------------------- Fig 7: perturbación vs modelo
fig, axs = plt.subplots(3, 3, figsize=(13, 8.4), sharex="col")
fig.subplots_adjust(top=0.85, bottom=0.10, left=0.07, right=0.99, hspace=0.12, wspace=0.22)
titulo(fig, "Paso 2: un océano que conduce explica buena parte de la perturbación medida",
       "Negro: perturbación medida (datos menos fondo). Naranja: campo de un dipolo inducido que se opone al campo de Júpiter, con la amplitud A\n"
       "ajustada. Lo que el modelo no explica (sobre todo en Bz) viene del plasma que choca con Europa, que este modelo todavía no incluye.")
for j, p in enumerate(PASOS):
    d, f = datos[p], fits[p]
    w = (d.r < 4).values
    mod = campo_dipolo(-0.5 * f["A"] * f["Bp"], d.x.values[w], d.y.values[w], d.z.values[w])
    for k, c in enumerate(["ΔBx", "ΔBy", "ΔBz"]):
        ax = axs[k, j]
        ax.plot(d["min"][w], f["dB"][w, k], color=INK, lw=0.9)
        ax.plot(d["min"][w], mod[k], color=C2, lw=2.2)
        ax.axhline(0, color="#c3c2b7", lw=0.8)
        if j == 0:
            ax.set_ylabel(f"{c} (nT)")
    axs[0, j].set_title(f"{p} · A ajustada = {f['A']:.2f}".replace(".", ","))
    axs[2, j].set_xlabel("Minutos desde la máxima aproximación")
axs[0, 0].legend(handles=[Line2D([], [], color=INK, lw=0.9, label="medido − fondo"),
                          Line2D([], [], color=C2, lw=2.2, label="dipolo inducido (modelo)")],
                 loc="lower left", frameon=False, fontsize=8.5)
nota(fig, "Ajuste por mínimos cuadrados con los datos a menos de 4 radios de Europa. Campo inductor = (Bx, By, 0) del fondo en la máxima aproximación (convención de Kivelson et al. 2000).")
fig.savefig(OUT + "fig7_perturbacion_vs_modelo.png", dpi=160); plt.close(fig)

# --------------------------------------------------- Fig 8: prueba de inversión + A
fig = plt.figure(figsize=(13, 5.6))
titulo(fig, "Paso 3: en E26 el campo de Júpiter se invierte, y el dipolo de Europa también",
       "Flechas en el plano ecuatorial (vista desde arriba). Azul: parte variable del campo de Júpiter. Naranja: dipolo de Europa ajustado libremente,\n"
       "sin imponer el modelo de inducción. En los tres pasos el dipolo apunta en contra del campo: eso es inducción, no un imán propio.")
gs = fig.add_gridspec(1, 4, left=0.04, right=0.99, top=0.76, bottom=0.12, wspace=0.35, width_ratios=[1, 1, 1, 1.35])
for j, p in enumerate(PASOS):
    ax = fig.add_subplot(gs[0, j]); ax.set_aspect("equal"); ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
    ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.add_patch(plt.Circle((0, 0), 0.45, color="#e9e2d4", ec="#c9bfae", zorder=1))
    f = fits[p]; Bp = f["Bp"][:2] / np.linalg.norm(f["Bp"][:2]); M = f["M"][:2] / np.linalg.norm(f["M"][:2])
    ax.annotate("", (1.25 * Bp[0], 1.25 * Bp[1]), (0, 0), arrowprops=dict(arrowstyle="-|>", color=C1, lw=3), zorder=3)
    ax.annotate("", (0.95 * M[0], 0.95 * M[1]), (0, 0), arrowprops=dict(arrowstyle="-|>", color=C2, lw=3), zorder=4)
    ang = np.degrees(np.arccos(np.clip(np.dot(M, -Bp), -1, 1)))
    ax.set_title(f"{p} · {FECHAS[p]}", loc="center")
    ax.text(0, -1.38, f"dipolo a {ang:.0f}° de lo esperado", ha="center", color=INK2, fontsize=9)
    ax.text(1.38, 0, "flujo →", ha="right", va="bottom", color=MUTED, fontsize=8)
    ax.text(0, 1.38, "↑ Júpiter", ha="center", va="top", color=MUTED, fontsize=8)
fig.legend(handles=[Line2D([], [], color=C1, lw=3, label="campo de Júpiter (parte variable)"),
                    Line2D([], [], color=C2, lw=3, label="dipolo de Europa (ajuste libre)")],
           loc="lower left", bbox_to_anchor=(0.04, 0.03), ncol=2, frameon=False, fontsize=9)
ax = fig.add_subplot(gs[0, 3])
for i, p in enumerate(PASOS):
    med, lo, hi = rob[p]
    ax.plot([lo, hi], [i, i], color=C2, lw=6, solid_capstyle="round", alpha=0.45)
    ax.plot(med, i, "o", color=C2, ms=9, mec="white", mew=1.5)
    ax.text(hi + 0.03, i, f"{med:.2f}".replace(".", ","), va="center", fontsize=9, color=INK)
ax.axvline(0.7, color=INK2, ls="--", lw=1); ax.text(0.69, -0.55, "mínimo exigido\npor Galileo (0,7)", fontsize=8, color=INK2, ha="right")
ax.axvline(1.0, color=C1, lw=1); ax.text(1.01, -0.55, "conductor\nperfecto", fontsize=8, color=C1)
ax.set_yticks(range(3)); ax.set_yticklabels(list(PASOS)); ax.set_xlim(0.5, 1.45); ax.set_ylim(-0.7, 2.6)
ax.invert_yaxis(); ax.set_xlabel("Amplitud A"); ax.set_title("A ajustada (27 variantes)")
nota(fig, "Barra: rango de A al variar el corte del fondo (3–4 RE), la ventana del ajuste (2,5–4 RE) y el grado del polinomio (2–4). Punto: mediana. "
     "A > 1 no es físico para un océano solo: indica la contribución del plasma.")
fig.savefig(OUT + "fig8_prueba_inversion.png", dpi=160); plt.close(fig)
print({p: tuple(round(v, 2) for v in rob[p]) for p in rob})
