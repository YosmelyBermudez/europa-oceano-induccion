# ¿Hay un océano bajo el hielo de Europa?
### Sondeo electromagnético de la luna Europa con Python

**Autora:** Yosmely Bermúdez, ingeniera geofísica (Universidad Simón Bolívar) y científica de datos
**Estado:** entrada 1 de la bitácora (14 nov 2026): modelo directo construido y validado
**Idioma:** español (habrá versión en inglés más adelante)

---

## 1. La pregunta, en una frase

¿Qué espesor y qué conductividad eléctrica puede tener el océano que hay bajo el hielo de Europa, según cómo responde la luna al campo magnético cambiante de Júpiter?

## 2. La idea: magnetotelúrica en una luna

En geofísica, la **magnetotelúrica (MT)** mide cómo responde el subsuelo a las variaciones naturales del campo magnético terrestre. Un campo que cambia induce corrientes eléctricas en las rocas que conducen. Esas corrientes generan un campo secundario que se mide en superficie. De la relación entre ambos campos se deduce la conductividad del subsuelo.

Europa vive exactamente esa situación:

| | Tierra (MT) | Europa |
|---|---|---|
| Fuente del campo variable | Sol + ionosfera | Júpiter: un imán inclinado ~10° que gira cada 9,925 h |
| Conductor | Roca con agua salada | Océano salado bajo el hielo |
| Sensor | Estación en superficie | Magnetómetro de una sonda (Galileo, Juno, Europa Clipper, JUICE) |
| Lo que se deduce | Conductividad del subsuelo | Espesor y conductividad del océano |

Como Europa orbita más lento de lo que gira Júpiter, el campo que siente oscila con el **periodo sinódico de 11,23 h**. Esa oscilación es la "señal de sondeo".

## 3. Antecedentes científicos (en qué se apoya)

- **Khurana et al. (1998)** interpretaron las perturbaciones magnéticas medidas por Galileo cerca de Europa y Calisto como campos inducidos por capas conductoras justo bajo la superficie, y argumentaron que lo más probable son océanos de agua salada.
- **Kivelson et al. (2000)** analizaron un sobrevuelo de enero de 2000 en el que el campo de Júpiter tenía la orientación contraria a la de los pasos anteriores. El dipolo de Europa también se invirtió. Eso distingue un campo inducido de uno permanente.
- **Zimmer, Khurana y Kivelson (2000)** acotaron la respuesta (A entre 0,7 y 1) y concluyeron que la conductividad del océano debe superar unos 0,058 S/m.
- **Schilling, Neubauer y Saur (2007)** añadieron la interacción con el plasma de Júpiter y obtuvieron conductividad × espesor ≥ 50 S/m·km.
- **Hand y Chyba (2007)** relacionaron esas cotas con la salinidad del océano.
- **Biersteker et al. (2023)** proponen una inversión bayesiana multifrecuencia pensada para los datos de Europa Clipper.
- **Connerney et al. (2022)** publicaron JRM33, el modelo del campo interno de Júpiter hecho con las 33 primeras órbitas de Juno. **Connerney et al. (2020)** publicaron Con2020, el modelo de la lámina de corriente (magnetodisco).

## 4. Estructura del proyecto

```
proyecto/
├── README.md                    ← este archivo
├── scripts/
│   ├── induccion.py             ← modelo de inducción (el corazón del proyecto)
│   ├── campo_inductor.py        ← campo de Júpiter en la órbita de Europa
│   ├── lectura_juno.py          ← lector de datos .sts de Juno
│   └── figuras_entrada01.py     ← figuras de la bitácora (entrada 1)
├── data/
│   ├── fgm_jno_l3_2016189pc_r60s_v01.sts   ← Juno FGM, 7 jul 2016, 60 s
│   └── B_europa_2per.npy        ← campo en Europa precalculado (2 periodos)
├── figs/                        ← 5 figuras de la entrada 1
└── carrusel/
    ├── estilo.py                ← estética "revista": renders de Júpiter y Europa
    ├── diapositivas.py          ← las 10 diapositivas del carrusel de LinkedIn
    ├── png/                     ← diapositivas 1080 × 1350 px
    └── Carrusel_Europa_oceano_LinkedIn.pdf
```

## 5. Qué hace cada código y por qué

### 5.1 `scripts/induccion.py`: la respuesta de Europa

**Qué hace.** Calcula la respuesta electromagnética Q = A·e^(iφ) de una esfera formada por capas concéntricas, cada una con su conductividad, ante un campo externo uniforme que oscila con periodo T.

**La física.** Sin corriente de desplazamiento, el campo magnético dentro de un conductor obedece la ecuación de difusión:

    ∇²B = μ₀ σ ∂B/∂t

Para un campo externo de grado 1 (uniforme), el campo dentro de cada capa se escribe con un escalar poloidal S(r):

- **En una capa conductora** (σ > 0): S = a·i₁(kr) + b·k₁(kr), donde i₁ y k₁ son funciones esféricas de Bessel modificadas y k = √(iωμ₀σ).
- **En una capa aislante** (σ = 0): S = a·r + b·r⁻².
- **En cada interfaz:** S y dS/dr son continuos.
- **Fuera de la luna:** S se descompone en una parte externa (el campo de Júpiter) y una interna (el dipolo inducido). Su cociente da Q.

**Convención (Zimmer et al. 2000):**

- **A = 1, φ = 0°:** conductor perfecto en la superficie. Expulsa todo el campo variable.
- **A → 0, φ → 90°:** mal conductor. Casi no responde.
- **Campo inducido fuera de la luna:** B_ind = −(Q/2)·(R/r)³·[3(B_p·r̂)r̂ − B_p].

**Funciones principales:**

- `respuesta(radios, sigmas, periodo_h)`: modelo general de N capas.
- `europa(h_oceano_km, sigma_oceano, hielo_km, periodo_h)`: atajo para el modelo de tres capas (interior aislante, océano, hielo).
- `esfera_uniforme_analitica(...)`: solución cerrada clásica, Q = 1 + 3/x² − 3·coth(x)/x con x = ka. Sirve de **validación independiente**.
- `respuesta_numerica(...)`: integra la ecuación diferencial con `scipy.integrate.solve_ivp` (DOP853). Es la **segunda validación independiente**.

### 5.2 `scripts/campo_inductor.py`: lo que siente Europa

**Qué hace.** Evalúa el campo de Júpiter en la órbita de Europa (9,38 RJ) a lo largo del tiempo y lo expresa en coordenadas EPhiO:

- x: dirección del flujo de plasma que corrota.
- y: hacia Júpiter.
- z: eje de giro de Júpiter.

**Modelos.** JRM33 hasta grado 18 (campo interno) más Con2020 (magnetodisco), ambos a través del paquete `JupiterMag`.

**Supuestos.** Órbita circular y ecuatorial. El periodo sinódico se calcula como 1/(1/9,925 h − 1/85,228 h) = 11,23 h.

**Resultado.**

- Componente hacia Júpiter: ±215 nT.
- Componente en la dirección del flujo: ±74 nT.
- Componente vertical: casi constante, alrededor de −396 nT.

### 5.3 `scripts/lectura_juno.py`: datos reales

**Qué hace.** Lee los archivos `.sts` del magnetómetro de Juno (PDS PPI, conjunto JNO-J-3-FGM-CAL-V1.0). Devuelve el tiempo, B (nT) y la posición (RJ) en el marco "planetocéntrico".

**Para qué se usó.** Para comprobar con datos reales que ese marco es System III y que el modelo JRM33 + Con2020 reproduce el campo medido (ver 6.2).

### 5.4 `scripts/figuras_entrada01.py`: figuras de la bitácora

Genera cinco figuras en `figs/`. Cada una lleva como título su conclusión, para que se entienda sin leer el texto.

1. `fig1_juno_marco_system3.png`: datos de Juno frente a los modelos.
2. `fig2_campo_primario_europa.png`: la señal de 11,23 h en Europa.
3. `fig3_validacion_modelo.png`: tres métodos que coinciden.
4. `fig4_espacio_oceanos.png`: mapa de amplitud A según espesor y conductividad.
5. `fig5_sobrevuelo_simulado.png`: señal esperada en un paso a 700 km.

### 5.5 `carrusel/`: divulgación

**`estilo.py`** reúne la estética "revista de astronomía":

- Fondo estelar.
- Render ortográfico de Júpiter coloreado con Br de JRM33 en la superficie, con iluminación y bandas.
- Render procedural de Europa (manchas y *lineae*).
- Corte de capas de Europa.
- Líneas de campo de un dipolo.

**`diapositivas.py`** genera las 10 diapositivas. Los gráficos se calculan con los mismos módulos científicos de `scripts/`, así que lo que se ve es lo que calcula el código.

## 6. Validación: números de control

La regla del proyecto es no usar nada sin antes reproducir un resultado conocido.

| # | Qué se comprueba | Resultado | Referencia | Estado |
|---|---|---|---|---|
| 1 | Código de capas frente a la solución analítica de la esfera uniforme | diferencia relativa < 10⁻¹⁴ | Solución cerrada clásica | ✓ |
| 2 | Código de capas frente a la integración numérica de la EDO | diferencia relativa < 10⁻¹¹ | Método independiente | ✓ |
| 3 | Límites físicos | A→1, φ→0° si σ→∞; A→0, φ→90° si σ→0 | Teoría | ✓ |
| 4 | Periodo de rotación del marco de los datos de Juno | 9,925 h | System III: 9,925 h | ✓ |
| 5 | Campo de Juno a 31–41 RJ frente a JRM33 + Con2020 | correlación ~0,96 (Bx, By); residuo rms ~5 nT | Connerney et al. (2020, 2022) | ✓ |
| 6 | Amplitud del campo inductor en Europa | 215 nT | ~200 nT (Biersteker et al. 2023) | ✓ |
| 7 | Conductividad mínima del océano para A ≥ 0,7 | **0,059 S/m** (con un océano de ~380 km) | **0,058 S/m** (Zimmer et al. 2000) | ✓ (2 %) |

### 6.1 Por qué los números 1 a 3 importan

Si tres métodos matemáticamente distintos dan lo mismo, el error no está en el código. Puede estar en los supuestos físicos, pero esos están documentados en la sección 7.

### 6.2 Cómo se obtuvo el número 4

La longitud de la sonda en el marco del archivo avanza −36,27°/h. Eso equivale a una vuelta cada 9,925 h, que es justo la rotación de System III. Juno casi no se mueve en ese tiempo, así que es el marco el que gira con Júpiter.

### 6.3 Cómo se obtuvo el número 7

Para cada conductividad σ se buscó el espesor del océano que da la amplitud máxima (`scipy.optimize.minimize_scalar`). Luego se buscó la σ para la que esa amplitud máxima vale 0,7 (`brentq`).

## 7. Supuestos y limitaciones (honestidad primero)

- **Simetría esférica y conductividad uniforme.** En la realidad la conductividad puede variar con la profundidad.
- **Sin plasma ni ionosfera.** La interacción con el plasma de Júpiter produce perturbaciones del mismo orden que la señal (Schilling et al. 2007). Habrá que modelarla o filtrarla antes de ajustar datos reales.
- **Solo el periodo de 11,23 h.** Falta el periodo orbital de 85,2 h y los armónicos.
- **Órbita de Europa circular y ecuatorial.**
- **El sobrevuelo de la figura 5 y de la diapositiva 8 es ilustrativo**, no un paso real de Galileo.
- **Los esquemas del carrusel no están a escala.** Las diapositivas 2 y 4 son ilustraciones conceptuales, no datos.

## 8. Cómo reproducirlo

```bash
pip install numpy scipy matplotlib pandas pillow JupiterMag

cd scripts
python figuras_entrada01.py           # figuras de la bitácora

cd ../carrusel
python diapositivas.py                # las 10 diapositivas en png/
python diapositivas.py 3 7            # solo algunas
```

Para regenerar el campo en Europa precalculado:

```bash
cd scripts
python -c "from campo_inductor import *; t=np.linspace(0,P_SIN*2,240,endpoint=False); B,_=campo_en_europa(t); np.save('../data/B_europa_2per.npy',np.vstack([t,B]))"
```

Para unir las diapositivas en el PDF del carrusel:

```bash
cd carrusel
python -c "from PIL import Image; ims=[Image.open(f'png/{i:02d}.png').convert('RGB') for i in range(1,11)]; ims[0].save('Carrusel_Europa_oceano_LinkedIn.pdf', save_all=True, append_images=ims[1:])"
```

## 9. Próximos pasos

1. Descargar del PDS PPI los datos del magnetómetro de Galileo en los sobrevuelos de Europa E4, E14 y E26.
2. Restar el campo de fondo y ajustar A y φ a los datos reales.
3. Repetir la prueba de Kivelson et al. (2000) con pasos de fases opuestas.
4. Inversión con incertidumbres (enfoque bayesiano al estilo de Biersteker et al. 2023).
5. Proyecto paralelo: profundidad de la fuente del campo de Júpiter con el espectro de Lowes. Connerney et al. (2022) reportan un radio de la dinamo de 0,81 RJ; mi estimación preliminar es ~0,80 RJ.

## 10. Referencias

1. Khurana, K. K., Kivelson, M. G., Stevenson, D. J., Schubert, G., Russell, C. T., Walker, R. J. y Polanskey, C. (1998). Induced magnetic fields as evidence for subsurface oceans in Europa and Callisto. *Nature*, 395, 777–780. https://doi.org/10.1038/27394
2. Kivelson, M. G., Khurana, K. K., Russell, C. T., Volwerk, M., Walker, R. J. y Zimmer, C. (2000). Galileo magnetometer measurements: A stronger case for a subsurface ocean at Europa. *Science*, 289.
3. Zimmer, C., Khurana, K. K. y Kivelson, M. G. (2000). Subsurface oceans on Europa and Callisto: Constraints from Galileo magnetometer observations. *Icarus*, 147, 329–347.
4. Schilling, N., Neubauer, F. M. y Saur, J. (2007). Time-varying interaction of Europa with the jovian magnetosphere: Constraints on the conductivity of Europa's subsurface ocean. *Icarus*, 192, 41–55.
5. Hand, K. P. y Chyba, C. F. (2007). Empirical constraints on the salinity of the europan ocean and implications for a thin ice shell. *Icarus*, 189, 424–438.
6. Biersteker, J. B., Weiss, B. P., Cochrane, C. J., Harris, C. D. K., Jia, X., Khurana, K. K., Liu, J., Murphy, N. y Raymond, C. A. (2023). Revealing the interior structure of icy moons with a Bayesian approach to magnetic induction measurements. *The Planetary Science Journal*, 4(4), 62. https://doi.org/10.3847/PSJ/acc331
7. Connerney, J. E. P., Timmins, S., Oliversen, R. J., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's prime mission. *Journal of Geophysical Research: Planets*, 127(2). https://doi.org/10.1029/2021JE007055
8. Connerney, J. E. P., et al. (2020). A Jovian magnetodisc model for the Juno era. *Journal of Geophysical Research: Space Physics*, 125.
9. Connerney, J. E. P. (2024). Juno MAG calibrated data, JNO-J-3-FGM-CAL-V1.0. NASA Planetary Data System, PPI node. https://doi.org/10.17189/1519711
10. James, M. JupiterMag: paquete de Python para los modelos de campo de Júpiter (PyPI).

**Agradecimiento a los datos:** este trabajo usa datos del Planetary Data System de la NASA (nodo PPI, UCLA), según su política de citación.
