"""Figuras de la entrada 1 de la bitácora."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import JupiterMag as jm
from lectura_juno import leer_sts
from induccion import *
from campo_inductor import campo_en_europa, P_SIN

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "lines.linewidth": 2})
OUT = "../figs/"


def titulo(fig, t, sub):
    fig.text(0.012, 0.985, t, fontsize=13.5, fontweight="bold", color=INK, va="top")
    fig.text(0.012, 0.935, sub, fontsize=10, color=INK2, va="top")


def nota(fig, txt):
    fig.text(0.012, 0.012, txt, fontsize=8, color=MUTED, va="bottom")


# ---------------------------------------------------------------- Fig 1: Juno
d = leer_sts("../data/fgm_jno_l3_2016189pc_r60s_v01.sts")
th = (d.t - d.t[0]).dt.total_seconds().values / 3600
jm.Internal.Config(Model="jrm33", CartesianIn=True, CartesianOut=True, Degree=18)
idx = np.arange(0, len(d), 10)
bi = np.array(jm.Internal.Field(d.x.values[idx], d.y.values[idx], d.z.values[idx]))
be = np.array(jm.Con2020.Field(d.x.values[idx], d.y.values[idx], d.z.values[idx]))
r = np.sqrt(d.x**2 + d.y**2 + d.z**2)
fig, axs = plt.subplots(3, 1, figsize=(10, 7.2), sharex=True)
fig.subplots_adjust(top=0.83, bottom=0.11, left=0.08, right=0.83, hspace=0.18)
titulo(fig, "Juno (7 jul 2016): el campo medido gira con Júpiter, y el marco «planetocéntrico» es System III",
       f"La posición de Juno en ese marco da una vuelta cada 9,925 h (periodo de System III). A {r.min():.0f}–{r.max():.0f} RJ "
       "el campo interno (JRM33) solo no basta:\nhace falta la lámina de corriente del magnetodisco (Con2020) para explicar los datos.")
for k, c in enumerate(["Bx", "By", "Bz"]):
    ax = axs[k]
    ax.plot(th, d[c], color=INK, lw=1.2)
    ax.plot(th[idx], bi[k], color=C1, ls="--")
    ax.plot(th[idx], bi[k] + be[k], color=C2)
    rr = np.sqrt(((d[["Bx", "By", "Bz"]].values[idx].T - (bi + be))[k]**2).mean())
    ax.set_ylabel(f"{c} (nT)")
    ax.text(1.01, 0.5, f"{c}\nresiduo rms\ncon disco: {rr:.1f} nT".replace(".", ","), transform=ax.transAxes, color=INK2, fontsize=9, va="center")
axs[-1].set_xlabel("Horas desde 2016-07-07 00:00 UTC (día 189)")
fig.legend(handles=[Line2D([], [], color=INK, lw=1.2, label="Juno FGM (medido)"),
                       Line2D([], [], color=C1, ls="--", label="JRM33 (solo interno)"),
                       Line2D([], [], color=C2, label="JRM33 + magnetodisco")],
              loc="upper left", bbox_to_anchor=(0.07, 0.885), ncol=3, frameon=False, fontsize=9.5)
nota(fig, "Datos: Juno FGM L3, PDS PPI (JNO-J-3-FGM-CAL-V1.0), 60 s. Modelos: JRM33 grado 18 (Connerney et al. 2022) y Con2020 (Connerney et al. 2020) vía JupiterMag.")
fig.savefig(OUT + "fig1_juno_marco_system3.png", dpi=160); plt.close(fig)

# ---------------------------------------------------- Fig 2: campo en Europa
t, Bx, By, Bz = np.load("../data/B_europa_2per.npy")
N = len(t); F = np.fft.rfft(np.vstack([Bx, By, Bz]), axis=1); amp = 2 * np.abs(F[:, 2]) / N
fig, axs = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
fig.subplots_adjust(top=0.84, bottom=0.12, left=0.09, right=0.84, hspace=0.12)
titulo(fig, f"Europa siente un campo que oscila ±{amp[1]:.0f} nT cada {P_SIN:.2f} h: esa es la señal que sondea su océano".replace("11.23", "11,23"),
       "Como el dipolo de Júpiter está inclinado ~10°, Europa pasa por encima y por debajo del ecuador magnético en cada periodo sinódico.\n"
       "Solo la parte que varía induce corrientes; la componente Bz es casi constante y no aporta.")
ax = axs[0]
ax.plot(t, Bx, color=C1); ax.plot(t, By, color=C2)
ax.axhline(0, color="#c3c2b7", lw=1)
ax.text(t[-1] + 0.3, Bx[-1], f"Bx (dirección del flujo)\namplitud {amp[0]:.0f} nT", color=INK2, fontsize=9, va="center")
ax.text(t[-1] + 0.3, By[-1], f"By (hacia Júpiter)\namplitud {amp[1]:.0f} nT", color=INK2, fontsize=9, va="center")
ax.annotate("", xy=(P_SIN + t[np.argmax(By)], 250), xytext=(t[np.argmax(By)], 250),
            arrowprops=dict(arrowstyle="<->", color=INK2, lw=1))
ax.text(t[np.argmax(By)] + P_SIN / 2, 262, f"periodo sinódico {P_SIN:.2f} h".replace(".", ","), ha="center", color=INK2, fontsize=9)
ax.set_ylim(-280, 300); ax.set_ylabel("nT")
ax = axs[1]
ax.plot(t, Bz, color=C3)
ax.text(t[-1] + 0.3, Bz[-1], f"Bz (eje de giro)\nmedia {Bz.mean():.0f} nT,\nvaría solo {amp[2]:.0f} nT", color=INK2, fontsize=9, va="center")
ax.set_ylabel("nT"); ax.set_xlabel("Horas (dos periodos sinódicos)")
nota(fig, "Modelo: JRM33 + Con2020 evaluados en la órbita de Europa (9,38 RJ, circular, ecuatorial). Ejes EPhiO. "
     "La literatura cita ~200 nT para el campo inductor sinódico.")
fig.savefig(OUT + "fig2_campo_primario_europa.png", dpi=160); plt.close(fig)

# --------------------------------------------------- Fig 3: validación modelo
sig = np.logspace(-3, 2, 60)
Qa = np.array([esfera_uniforme_analitica(R_EUROPA, s, P_SIN) for s in sig])
Qb = np.array([respuesta([R_EUROPA], [s], P_SIN) for s in sig])
sn = sig[::6]
Qn = np.array([respuesta_numerica([R_EUROPA], [s], P_SIN) for s in sn])
# también para un océano de 100 km bajo 25 km de hielo
r_top = R_EUROPA - 25e3
Qo = np.array([europa(100, s, 25, P_SIN) for s in sig])
Qon = np.array([respuesta_numerica([r_top - 100e3, r_top, R_EUROPA], [0, s, 0], P_SIN) for s in sn])
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.9), gridspec_kw={"width_ratios": [1, 1, 0.9]})
fig.subplots_adjust(top=0.78, bottom=0.14, left=0.06, right=0.98, wspace=0.32)
err = max(np.max(np.abs(Qa - Qb) / np.abs(Qa)), np.max(np.abs(Qn - Qb[::6]) / np.abs(Qn)),
          np.max(np.abs(Qon - Qo[::6]) / np.abs(Qon)))
titulo(fig, f"El código de inducción coincide con dos soluciones independientes (diferencia máxima {err:.0e})".replace("e-", "×10^-").replace("×10^-12", "×10⁻¹²").replace("×10^-11", "×10⁻¹¹"),
       "Línea: solución cerrada clásica (esfera uniforme). Puntos: mi código de capas con funciones de Bessel. Cruces: integración numérica de la EDO.\n"
       "Límites físicos correctos: un mal conductor no responde (A→0, desfase→90°) y un conductor perfecto expulsa el campo (A→1, desfase→0°).")
for ax, f, lab in [(axs[0], np.abs, "Amplitud A"), (axs[1], lambda q: np.degrees(np.angle(q)), "Desfase φ (°)")]:
    ax.plot(sig, f(Qa), color=C1, label="Esfera uniforme: analítica")
    ax.plot(sig[::3], f(Qb[::3]), "o", ms=5, mfc="none", mec=C1)
    ax.plot(sn, f(Qn), "x", ms=8, color=INK)
    ax.plot(sig, f(Qo), color=C2)
    ax.plot(sig[::3], f(Qo[::3]), "o", ms=5, mfc="none", mec=C2)
    ax.plot(sn, f(Qon), "x", ms=8, color=INK)
    ax.set_xscale("log"); ax.set_xlabel("Conductividad σ (S/m)"); ax.set_title(lab)
axs[0].text(1.1e-3, 0.42, "esfera\nuniforme", color=C1, fontsize=9)
axs[0].text(0.6, 0.35, "océano de 100 km\nbajo 25 km de hielo", color=C2, fontsize=9)
axs[0].legend(handles=[Line2D([], [], color=INK2, label="solución analítica"),
                       Line2D([], [], ls="", marker="o", mfc="none", mec=INK2, label="mi código (capas)"),
                       Line2D([], [], ls="", marker="x", color=INK, label="EDO numérica")], loc="upper left", frameon=False, fontsize=8.5)
ax = axs[2]
ax.semilogy(sig, np.abs(Qa - Qb) / np.abs(Qa) + 1e-17, color=C1, label="capas vs analítica")
ax.semilogy(sn, np.abs(Qn - Qb[::6]) / np.abs(Qn) + 1e-17, "x-", color=INK, lw=1, label="capas vs EDO (esfera)")
ax.semilogy(sn, np.abs(Qon - Qo[::6]) / np.abs(Qon) + 1e-17, "x-", color=C2, lw=1, label="capas vs EDO (océano)")
ax.set_xscale("log"); ax.set_title("Diferencia relativa |ΔQ|/|Q|"); ax.set_xlabel("Conductividad σ (S/m)")
ax.legend(frameon=False, fontsize=8.5, loc="upper left"); ax.set_ylim(1e-17, 1e-6)
nota(fig, f"Periodo del forzamiento: {P_SIN:.2f} h.".replace("11.23","11,23") + f" Radio de Europa: 1560,8 km. Q = A·e^(iφ), con A = 1 para un conductor perfecto en la superficie (convención de Zimmer et al. 2000).")
fig.savefig(OUT + "fig3_validacion_modelo.png", dpi=160); plt.close(fig)

# ------------------------------------------------- Fig 4: espacio de océanos
hs = np.linspace(5, 600, 90); ss = np.logspace(-3, np.log10(30), 90)
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True)
fig.subplots_adjust(top=0.80, bottom=0.12, left=0.07, right=0.90, wspace=0.08)
titulo(fig, "Para dar A ≥ 0,7 (Galileo), el océano necesita σ ≥ 0,059 S/m. Zimmer et al. (2000) obtuvieron 0,058 S/m",
       "Cada punto es un océano posible (espesor y conductividad). Color: amplitud de la respuesta a 11,23 h. Línea gruesa: A = 0,7, el mínimo que exigen\n"
       "los datos de Galileo según Zimmer et al. (2000); solo los océanos a la derecha de ella son compatibles. Línea discontinua: desfase en grados.")
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("azul", ["#f0efec", "#cde2fb", "#86b6ef", "#2a78d6", "#184f95", "#0d366b"])
for ax, hielo in zip(axs, [0, 25]):
    Q = np.array([[europa(h, s, hielo, P_SIN) for s in ss] for h in hs])
    A, ph = np.abs(Q), np.degrees(np.angle(Q))
    cf = ax.contourf(ss, hs, A, levels=np.linspace(0, 1, 11), cmap=cmap)
    cs = ax.contour(ss, hs, A, levels=[0.7], colors=C2, linewidths=3)
    cp = ax.contour(ss, hs, ph, levels=[5, 10, 20, 45], colors=INK2, linewidths=0.8, linestyles="--")
    ax.clabel(cp, fmt="%d°", fontsize=8)
    ax.axvspan(2.75, 5, color=C3, alpha=0.18, lw=0)
    ax.text(3.6, 590, "agua de mar\nterrestre", color="white", fontsize=8.5, ha="center", va="top")
    ax.set_xscale("log"); ax.set_xlabel("Conductividad del océano σ (S/m)")
    ax.set_title(f"Corteza de hielo de {hielo} km" + (" (caso de Zimmer)" if hielo == 0 else ""))
    ax.grid(False)
    if hielo == 0:
        ax.plot(0.0591, 380, "o", ms=9, color=C2, mec="white", mew=1.5)
        ax.annotate("mínimo: σ = 0,059 S/m\ncon océano de ~380 km", (0.0591, 380), (0.0025, 470), fontsize=9, color=INK,
                    arrowprops=dict(arrowstyle="->", color=INK2))
        ax.text(0.4, 60, "A = 0,7", color=C2, fontsize=10, fontweight="bold")
axs[0].set_ylabel("Espesor del océano h (km)")
cb = fig.colorbar(cf, ax=axs, fraction=0.025, pad=0.02); cb.set_label("Amplitud A")
nota(fig, "Modelo de 3 capas: interior aislante, océano conductor uniforme y corteza de hielo aislante. Sin ionosfera ni interacción con el plasma.")
fig.savefig(OUT + "fig4_espacio_oceanos.png", dpi=160); plt.close(fig)

# ------------------------------------------------------ Fig 5: sobrevuelo
Bhat = F[:, 2] * 2 / N                               # fasor del campo primario (nT)
i0 = np.argmax(By); w = 2 * np.pi / P_SIN
casos = [("Sin océano", None, INK2, ":"),
         ("Océano mínimo\n  σ = 0,06 S/m, h = 380 km, sin hielo", europa(380, 0.06, 0, P_SIN), C1, "-"),
         ("Agua de mar\n  σ = 3 S/m, h = 100 km, 25 km de hielo", europa(100, 3.0, 25, P_SIN), C2, "-")]
R = R_EUROPA / 1e3; alt = 700.0; v = 6.0                # km, km/s
tm = np.linspace(-15, 15, 400); s = tm * 60 * v
pos = np.stack([s, -(R + alt) * np.ones_like(s), 0 * s])  # paso por el lado opuesto a Júpiter
rn = np.linalg.norm(pos, axis=0); rh = pos / rn
fig, axs = plt.subplots(3, 1, figsize=(10, 7.4), sharex=True)
fig.subplots_adjust(top=0.84, bottom=0.13, left=0.09, right=0.68, hspace=0.15)
titulo(fig, "Un sobrevuelo a 700 km distinguiría un océano «mínimo» de uno tipo agua de mar",
       "Campo inducido esperado (con el campo de Júpiter ya restado) en un paso rectilíneo simulado, cuando By es máximo.\n"
       "Los dos océanos difieren en 10–25 nT en Bx y By, muy por encima del ruido del magnetómetro. Bz casi no cambia.")
for lab, Q, col, ls in casos:
    if Q is None:
        Bi = np.zeros_like(pos)
    else:
        M = -0.5 * np.real(Q * Bhat * np.exp(1j * w * t[i0]))
        M[2] = -0.5 * np.real(Q * Bhat[2] * np.exp(1j * w * t[i0]))
        Bi = (R / rn)**3 * (3 * (M[:, None] * rh).sum(0) * rh - M[:, None])
        lab = lab + f"\n  → A = {abs(Q):.2f}, φ = {np.degrees(np.angle(Q)):.0f}°".replace(".", ",")
    for k in range(3):
        axs[k].plot(tm, Bi[k], color=col, ls=ls, label=lab if k == 0 else None)
for k, c in enumerate(["ΔBx", "ΔBy", "ΔBz"]):
    axs[k].set_ylabel(f"{c} (nT)")
axs[-1].set_xlabel("Minutos desde la máxima aproximación")
axs[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8.5)
nota(fig, "Geometría ilustrativa (no es un sobrevuelo real): paso a 6 km/s en la dirección del flujo, a 700 km de altura (como E4).\n"
     "Solo dipolo inducido: no incluye las perturbaciones del plasma, que en los datos reales son del mismo orden.")
fig.savefig(OUT + "fig5_sobrevuelo_simulado.png", dpi=160); plt.close(fig)
print("ok")
