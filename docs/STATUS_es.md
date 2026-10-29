# AFA-12 v3 — Estado, revisión del diseño original y trabajo previo

> **Nota:** AFA-12 fue el nombre de trabajo del proyecto. Se publica como **Gosset-KEM** (especificación 1.0). Los juegos AFA12-… de este documento corresponden a Gosset-KEM-…; las etiquetas internas de la versión publicada son distintas, por eso sus vectores de prueba difieren de los citados acá.

**Fecha:** 2026-10-24

## 1. Reordenamiento del diseño

| Ley / componente (Flujo V6) | v1 | v3 |
|---|---|---|
| Ley 1 — Materia y Fuerza (Binary Embedding) | m_i → m_i·⌈q/2⌋ (1 bit por coeficiente) | **Codificación E8: 8 bits cada 8 coeficientes** |
| Ley 3 — Masa de Higgs (ruido) | Gaussiana sobre E8^32 | CBD estándar (óptimo, ver v2) |
| Ley 6 — Fermión Fantasma (secreto) | Sobre E8 | CBD estándar |
| Ley 2, 4, 5, FO, DEM | sin cambios | ver puntos pendientes |

Resultado exacto (`v3/e8_message_coding.py`, validado contra la especificación de Kyber): con codificación E8 la tasa de fallo de ML-KEM-512 baja de 2^-139 a 2^-241. Con u comprimido a 9 bits: ciphertext de 704 bytes (−64 bytes) con fallo 2^-142.

## 2. Puntos del diseño original pendientes de revisar

| # | Dónde | Observación | Riesgo |
|---|---|---|---|
| 1 | Paper: n=256, A ∈ Z_q^{256×256} (LWE plano) | Con el estimador validado da unos 45 bits clásicos, muy lejos de NIST nivel 1 (~118). Hace falta estructura de módulo k≥2, como ML-KEM. | Alto |
| 2 | Flujo V6, Ley 2: "r ∈ Z_q^256 vía SHAKE-256" | Si r fuera uniforme mod q, el término eᵀ·r sería enorme y el descifrado fallaría siempre. Tiene que ser chico (CBD), como en Kyber. Verificar si es solo notación. | Alto si es literal |
| 3 | Flujo V6, capa DEM | La misma clave K se usa para el keystream (SHAKE) y para HMAC. Conviene derivar dos claves independientes con un KDF. | Medio |
| 4 | Paper sección 7 | El texto extraído perdió las fórmulas porque son imágenes (120). Conviene reescribirlas en LaTeX o con el editor de ecuaciones. | Documental |
| 5 | Nombres de las leyes (Higgs, Fermión, etc.) | Sirven como metáfora interna; para publicar hay que usar la terminología criptográfica estándar. | Bajo |
| 6 | Paper sección 13 (trabajo futuro: Leech, D4, Barnes-Wall, anisotrópico) | Como **ruido**, el resultado de la v2 ya los descarta (cualquier sublattice pierde log₂D bits). Como **código del mensaje** siguen abiertos. | — |

## 3. Trabajo previo sobre codificación con lattices en Kyber

| Trabajo | Qué hace | Qué deja abierto |
|---|---|---|
| Saliba, Luzzi, Ling — E8 (Gosset) para FrodoKEM, arXiv:2110.01740 | Códigos de 8 dimensiones; requiere cambiar parámetros | — |
| Liu & Sakzad — Lattice Codes for CRYSTALS-Kyber, arXiv:2308.13981 | BW16 y Leech24 en Kyber: 20–33% menos de ciphertext; decodificadores argumentados como de tiempo constante | No hay mediciones de side-channel; el DFR viene de una aproximación por TCL |
| eprint 2024/1243 — códigos 2D para Kyber | Códigos 2D sin cambiar parámetros; implementación isócrona (404 ciclos de decodificación) | Ganancia chica (Nivel 5: −6% de ciphertext) |
| **Liu & Sakzad — Maximum Toroidal Distance Code, arXiv:2601.08452 (enero 2026)** | **E8 para Kyber (código GTD). Kyber1024: DFR 2^-286 (du=11, dv=5)** | **El decodificador E8 no está implementado ni medido; los resultados numéricos son solo para Kyber1024; DFR por cota de unión, sin simulación** |

**Conclusión:** la codificación E8 en Kyber ya está publicada (enero 2026). El aporte posible de AFA-12 está en los huecos:
1. Un decodificador E8 de tiempo constante integrado en ML-KEM, con benchmarks reales y pruebas de side-channel.
2. Resultados para ML-KEM-512 y 768, que el paper de 2026 no cubre en detalle.
3. Validar el DFR más allá de la cota de unión y del supuesto de independencia.
4. ~~Reconciliar nuestro cálculo (2^-335) con el de ellos (2^-286).~~ **Resuelto, ver sección 4.**

## 4. Reconciliación con Liu & Sakzad 2026 (`v3/reconcile.py`)

La diferencia es de modelo de ruido, no de error. Con la misma cota de unión:

| Caso | var(w) | Estándar exacto | Estándar gaussiano | E8 exacto | E8 gaussiano | Liu & Sakzad |
|---|---|---|---|---|---|---|
| ML-KEM-512 (10,4) | 6238 | −139.1 | −76.0 | −240.7 | −152.2 | — |
| ML-KEM-768 (10,4) | 5854 | −165.3 | −81.3 | −279.4 | −162.7 | — |
| ML-KEM-1024 (11,5) | 3345 | −175.2 | −145.8 | −335.0 | −290.9 | −286 |
| ML-KEM-1024 (10,5) | 3897 | −143.3 | −124.5 | −275.7 | −248.6 | −239 |

(valores en log₂ de la tasa de fallo)

- Con la aproximación gaussiana (teorema central del límite) reproducimos sus cifras con unos 5–10 bits de diferencia, atribuibles a detalles de la cota o del decodificador.
- El cálculo exacto (convolución) es el mismo método que usa la especificación oficial de Kyber, y reproduce los valores publicados (−139/−165/−175). La aproximación gaussiana es muy pesimista cuando dv=4, porque el error de compresión es acotado y no gaussiano.
- **Para comparar contra las cifras de NIST, el método correcto es el exacto.** Ellos subestiman la ganancia de E8.

**Supuesto que quedaba abierto:** tanto nuestro cálculo como el de ellos tratan como independientes los ruidos de los 8 coeficientes de un bloque, aunque en un cifrado real comparten s y e. Se validó en la sección 5.

## 5. Validación con cifrados reales (`v3/correlation_check.py`, `v3/tail_check.py`)

Se usó un PKE tipo ML-KEM completo: anillo Z_q[X]/(X²⁵⁶+1), rango de módulo k, CBD, compresión du=10 y dv=4, con η inflado para que los fallos sean observables. Hubo claves nuevas en cada lote (entre 400 y 6000 claves) y entre 1,28 y 19,2 millones de bloques de 8 coeficientes por configuración. Se compara la tasa de fallo por bloque medida (real) contra coeficientes independientes con la misma distribución marginal exacta (iid) y contra la cota de unión.

| k | η | Codificación | Bloques | Real | iid | Cota de unión | Real/iid |
|---|---|---|---|---|---|---|---|
| 2 | 12 | estándar | 19,2 M | 4,995e-4 | 4,984e-4 | 4,939e-4 | 1,00 |
| 2 | 12 | **E8** | 19,2 M | 1,875e-6 (36 ev.) | 1,550e-6 (31 ev.) | 1,869e-6 | 1,21 (n.s.) |
| 2 | 13 | estándar | 12,8 M | 1,559e-3 | 1,549e-3 | 1,551e-3 | 1,01 |
| 2 | 13 | **E8** | 12,8 M | 1,711e-5 (219 ev.) | 1,600e-5 (208 ev.) | 1,705e-5 | 1,07 (n.s.) |
| 2 | 14 | estándar | 1,28 M | 4,100e-3 | 4,002e-3 | 3,951e-3 | 1,02 |
| 2 | 14 | **E8** | 1,28 M | 1,016e-4 | 1,125e-4 | 1,037e-4 | 0,90 |
| 2 | 16 | estándar | 1,28 M | 1,670e-2 | 1,615e-2 | 1,639e-2 | 1,03 |
| 2 | 16 | **E8** | 1,28 M | 1,445e-3 | 1,352e-3 | 1,618e-3 | 1,07 |
| 2 | 18 | estándar | 1,28 M | 4,542e-2 | 4,516e-2 | 4,520e-2 | 1,01 |
| 2 | 18 | **E8** | 1,28 M | 8,764e-3 | 8,235e-3 | 1,155e-2 | 1,06 |
| 3 | 16 | estándar | 1,28 M | 8,654e-2 | 8,596e-2 | 8,595e-2 | 1,01 |
| 3 | 16 | **E8** | 1,28 M | 2,579e-2 | 2,523e-2 | 4,063e-2 | 1,02 |

(n.s. = diferencia no significativa con los eventos observados)

**Conclusiones:**
- **El supuesto de independencia se sostiene** en todo el rango medido, para ambas codificaciones. Hay un posible exceso de alrededor del 5% en E8 con ruido moderado, que equivale a menos de 0,1 bit. Al bajar el ruido no se detecta que crezca: en η=12 y 13 la diferencia es estadísticamente nula.
- **En la cola, la cota de unión coincide casi exactamente con la tasa real medida:** E8 en η=13 predice 1,705e-5 y se mide 1,711e-5; en η=12, 1,869e-6 contra 1,875e-6. Las cifras criptográficas de la sección 4 descansan sobre una cota que, en la zona medible, es prácticamente exacta.
- **Límite:** no se puede simular hasta 2⁻²⁴⁰. La extrapolación desde 10⁻⁶ es el mismo supuesto que acepta la especificación de Kyber para sus propias cifras.

## 6. v4 — Decodificador E8 de tiempo constante en C (`v4/`)

**Implementación** (`e8codec.c`, unas 150 líneas):
- **Codificación:** 8 bits por bloque de 8 coeficientes, como punto de (q/2)·E8 mod q·ℤ⁸. Los valores posibles son {0, 832, 1665, 2497}.
- **Decodificación:** algoritmo de Conway–Sloane (dos candidatos, D8 y D8+½, con corrección de paridad), todo en aritmética entera.
- **Sin dependencias del secreto:** no hay ramas que dependan de datos secretos ni accesos a memoria indexados por ellos. El redondeo se hace por multiplicación y desplazamiento, sin instrucción de división, y las selecciones con máscaras.
- **Ensamblador revisado** (gcc 13, -O2 y -O3): cero `div`/`idiv`. Los únicos saltos condicionales son contadores de bucle con límites constantes y el chequeo del stack protector.

**Corrección** (`test_correctness.py`):

| Prueba | Resultado |
|---|---|
| Redondeo por multiplicación y desplazamiento, todo el rango (106.529 entradas) | 0 diferencias |
| Ida y vuelta de los 256 mensajes | OK |
| Contra el decodificador Python, 2 M bloques con entrada uniforme | 391 diferencias, **todas empates exactos** de distancia (ambas respuestas válidas) |
| Fallo con ruido ML-KEM real (k=2, η=13, 4,8 M bloques) | 2,0e-5 (96 eventos), compatible con 1,71e-5 de la v3 |

**Benchmark** (Xeon 2,1 GHz, mediana de ciclos por mensaje de 256 coeficientes):

| Función | Ciclos |
|---|---|
| Decodificación estándar ML-KEM (`poly_tomsg`) | 1.574 |
| **Decodificación E8 de tiempo constante** | **5.420 (3,4×)** |
| Codificación E8 | 638 |

El costo extra es de unos 3.800 ciclos por desencapsulación. Su peso relativo frente a una desencapsulación ML-KEM completa queda por medir con una implementación entera.

**Corrección de la línea base y versión optimizada** (`e8fast.c`, `bench_fast.c`, gcc -O3 -mavx2):

La línea base anterior (1.574 ciclos) era una versión lenta propia de `poly_tomsg`, con el operador `%`. La función de referencia de Kyber cuesta 274 ciclos, así que el sobrecosto real de la v4.0 era de 16,9×, no de 3,4×. La versión optimizada aprovecha que |Y| < U: el redondeo se reduce a 2 comparaciones, sin multiplicaciones ni divisiones, y toda la aritmética es int32.

| Decodificador | Ciclos | vs `poly_tomsg` |
|---|---|---|
| `poly_tomsg` de referencia ML-KEM | 274 | 1,0× |
| E8 v4.0 (primera versión) | 4.638 | 16,9× |
| **E8 rápido (escalar int32)** | **1.030** | **3,8×** |
| E8 AVX2 (un bloque por registro) | 1.266 | 4,6× |

- **Corrección:** la versión rápida y la AVX2 dan exactamente la misma salida que el decodificador verificado en 6,4 millones de bloques (0 diferencias).
- **Tiempo constante** (5 M mediciones): rápido |t| = 1,30, AVX2 |t| = 1,34, referencia |t| = 1,24, sin fuga detectada en ninguno. El control con fugas da |t| ≈ 23.300.
- **Costo extra real:** unos 760 ciclos por desencapsulación. La desencapsulación completa de ML-KEM cuesta decenas de miles de ciclos, así que el peso relativo debería ser de pocos por ciento. Falta confirmarlo con una implementación completa.
- La versión AVX2 por bloque resultó más lenta que la escalar, porque las operaciones horizontales (máximo, paridad, suma) dominan. Una versión que procese 8 bloques en paralelo (una coordenada por registro) evitaría ese costo; queda como optimización pendiente.
- **Endurecimiento posterior** (lección del decodificador Leech): se agregó la misma barrera de valor a `e8fast.c`. La versión escalar pierde la vectorización automática y pasa a 4.418 ciclos; **la versión AVX2 (intrínsecos) no cambia: 1.346 ciclos (5,1×), |t| = 2,27**, y es la recomendada. La escalar sin barrera (1.030 ciclos) estaba limpia con gcc 13, pero no hay garantía con otros compiladores. Corrección revalidada: 0 diferencias en 6,4 M bloques.

**Prueba de tiempo constante** (`bench_ct.c`, metodología dudect: test t de Welch, clase fija contra aleatoria, recorte por percentiles, entradas preparadas antes de medir):

| Función | Mediciones | máx \|t\| | Veredicto |
|---|---|---|---|
| e8_decode | 10 M (3 corridas) | 1,36 / 2,01 / 1,56 | sin fuga detectada |
| std_decode | 10 M | 0,93 – 1,69 | sin fuga detectada |
| e8_decode_LEAKY (control positivo con ramas intencionales) | 10 M | ~26.000 | **fuga detectada** |

Lección metodológica: una primera versión del arnés preparaba cada entrada de forma distinta según la clase, justo antes de medir. Eso producía falsos positivos (|t| hasta 8) incluso en la función estándar. Se corrigió preparando cada lote de entradas por adelantado.

**Límites:** se midió en una VM cloud con un solo compilador y una sola CPU, y solo el canal de tiempo. No se evaluaron canales de potencia ni electromagnéticos, ni otras plataformas (ARM, microcontroladores).

## 7. Revisión del documento "AFA-12 v2.0" (22-jul-2026, generado con otra herramienta)

Ese documento propuso subir los parámetros a n=768 y σ=32.5 para alcanzar "128 bits". Se verificó con el estimador validado y el modelo de ruido:

| Afirmación del v2.0 | Verificación | Veredicto |
|---|---|---|
| β = n·log(q/σ)/log(q) ≈ 438 → 128 bits | Esa fórmula no corresponde a ningún estimador estándar. Con el estimador validado (que reproduce Kyber) y n=768, σ=32.5 (desvío 12,97), da β≈1488 (~434 bits) | Fórmula incorrecta |
| Con n=768 y σ=32.5 el esquema funciona | sd(w) = 6.588 frente a q/4 = 832: falla el 90% de los coeficientes y prácticamente el 100% de los mensajes. Aun con s y r chicos (CBD2), sd(w) = 508 y falla el 10% de los coeficientes | **El esquema no descifra** |
| "D_{E8,σ} mod q es ε-cercana a la uniforme (Δ ≤ 2⁻¹²⁹)" | Distancia de variación total medida: 0,974. Además, si el ruido fuera uniforme mod q el descifrado sería imposible. El parámetro de suavizado garantiza uniformidad módulo E8, no módulo q | Mal aplicado |
| "Reducción worst-case parcial vía GLWE" | Se apoya en la afirmación anterior | No válida |
| "Todas las normas de Gram-Schmidt valen 1" | Normas reales: [2, 1, 1, 1, 1, 1, 1, 0,5] | Incorrecto |
| ρ_max ≈ 0,925 entre coordenadas adyacentes | Medido con el sampler correcto: \|ρ\| < 0,02 | Incorrecto |
| Varianza de CBD η=2 = 0,5 | Es 1,0 | Incorrecto |
| Overhead del muestreo ~24× | Plausible en Python; irrelevante en la v3 (el ruido vuelve a ser CBD) | — |

**Lo rescatable:** la intuición de "llevar n a 768" es correcta, si se hace con ruido chico y con estructura de módulo. Con n_total=768 y CBD2, el estimador validado da β=624 (182 bits). Eso es exactamente ML-KEM-768 (k=3 × 256), que además es rápido gracias a la NTT sobre polinomios, en vez de una matriz 768×768 sin estructura. **El diseño v3 (ML-KEM + codificación E8) ya incorpora ese ajuste.**

### 7.1 ¿Existe otro σ que haga funcionar el diseño v2.0? (`v6/sigma_sweep.py`, `v6/sigma_max.py`)

Diseño del v2.0: LWE plano, n=768, A uniforme 768×768, y s, e, r, e₁, e₂ gaussianas discretas de desvío σ, sin compresión. Se barrió σ calculando la seguridad con el estimador validado y la tasa de fallo exacta por convolución.

| q | σ (desvío) | Bits, ruido estándar | Bits, ruido E8 | log₂ fallo, código estándar | log₂ fallo, código E8 |
|---|---|---|---|---|---|
| 3329 | 0,8 | 171 | 156 | −632 | < −999 |
| 3329 | 1,0 | 182 | 166 | −287 | −573 |
| 3329 | 1,2 | 192 | 174 | −144 | −288 |
| 3329 | 1,5 | 204 | 185 | −58 | −118 |
| 3329 | 2,0 | 222 | 201 | −15 | −32 |
| 3329 | 12,97 (v2.0) | 434 | 376 | +8 (falla siempre) | +12 |
| 2¹⁵ | 2,8 | 177 | 162 | −432 | −865 |
| 2¹⁵ | 4,0 | 195 | 178 | −113 | −228 |
| 2¹⁵ | 6,0 | 219 | 199 | −19 | −40 |

**σ máximo con fallo ≤ 2⁻¹²⁸** (q=3329, n=768, búsqueda binaria exacta):
- **Código estándar:** σ_max = 1,235 → 193 bits (ruido estándar) / 175 bits (ruido E8).
- **Código E8:** σ_max = 1,470 → **203 bits** (ruido estándar) / 184 bits (ruido E8).

**Conclusiones:**
- **El σ correcto existe y está entre 1 y 1,5,** unas 10 veces menos que 32.5 en desvío. Es la misma zona que usa Kyber (CBD η=2 → σ=1, η=3 → σ≈1,22).
- **Otra forma de cobrar la ganancia de E8:** en vez de achicar el ciphertext, se puede subir σ de 1,235 a 1,47 y ganar unos **10 bits de seguridad** con el mismo fallo.
- **El ruido E8 sigue perdiendo** entre 15 y 19 bits frente al ruido estándar con el mismo σ, en todos los casos.
- **El requisito de suavizado del v2.0** (σ_ρ > 7,63, o sea desvío > ~3,0) es incompatible con la corrección cuando q=3329. Solo cabe con q=2¹⁵ (σ hasta ~4), que es el terreno de FrodoKEM (n=640, q=2¹⁵, σ=2,8). Saliba, Luzzi y Ling ya aplicaron E8 como código del mensaje en FrodoKEM.
- **Velocidad:** el LWE plano del v2.0 necesita multiplicar matrices de 768×768 (~590.000 productos por operación), mientras que ML-KEM usa la NTT sobre polinomios. Por eso FrodoKEM es mucho más lento que ML-KEM. El v2.0 va en contra del objetivo de velocidad.

### 7.2 Resto del v2.0 que no se había probado

| Punto | Estado |
|---|---|
| Muestreo anisotrópico (trabajo futuro 3) | No probado experimentalmente. Por el mismo argumento de la v2, con varianza total fija (lo que fija la corrección), la entropía y el volumen efectivo del ruido son máximos cuando el ruido es isotrópico (AM-GM sobre los autovalores). No se espera ganancia. |
| Híbrido Module-LWE + ruido Λ (trabajo futuro 5) | Cubierto por el teorema de la v2: el ruido en un sublattice pierde log₂D bits, sea cual sea la estructura de A. |
| Cálculo exacto de η_ε con la función theta | Irrelevante: el requisito de suavizado no se usa en el diseño viable. |
| Arora-Ge | Con ruido gaussiano o CBD de varianza ≥ 1 no es práctico (igual que en Kyber). Sin cambios. |
| Capa DEM | Pendiente: derivar claves separadas para el keystream y para el HMAC (hoy se reusa K). |
| Fase 5 del v2.0 (SIMD para SampleE8) | Ya no aplica: el ruido es CBD y el componente SIMD es el decodificador E8 (AVX2, 1.346 ciclos). |

## 8. v5 — Leech (Λ₂₄) contra E8 contra el estándar (`v5/`)

**Por qué Leech no sirve a 1 bit por coeficiente:** con 1 bit por coeficiente, un lattice L con qℤ^ℓ ⊂ L e índice 2^ℓ es (q/2) veces un lattice unimodular que contiene 2ℤ^ℓ. Eso equivale a una Construcción A de un código binario autodual, cuya norma mínima es como máximo 2: E8 ya alcanza el óptimo. La ventaja de Leech aparece a **1,5 bits por coeficiente**: L = (q/8)·Λ₂₄ (coordenadas enteras, norma mínima² 32, contiene 8ℤ²⁴) tiene 2³⁶ puntos por cada 24 coeficientes y **la misma distancia mínima que E8 (0,707·q)**. Un mensaje de 256 bits entra en 8 bloques, es decir 192 coeficientes de v; los 64 restantes no se transmiten. A cambio, tiene 196.560 vecinos mínimos contra 240.

**Cálculo exacto** (`leech_eval.py`, mismo método validado contra la especificación de Kyber; cota de unión por tipo de vector mínimo: 1104 (±4², 0²²), 97152 (±2⁸, 0¹⁶), 98304 (∓3, ±1²³)):

| Nivel | Fallo estándar | Fallo E8 | Fallo Leech (mismos du, dv) |
|---|---|---|---|
| ML-KEM-512 | 2⁻¹³⁹·¹ | 2⁻²⁴⁰·⁷ | 2⁻¹⁹⁹·⁹ |
| ML-KEM-768 | 2⁻¹⁶⁵·³ | 2⁻²⁷⁹·⁴ | 2⁻²²⁸·² |
| ML-KEM-1024 | 2⁻¹⁷⁵·² | 2⁻³³⁵·⁰ | 2⁻³¹³·¹ |

**Ciphertext mínimo que cumple el objetivo de fallo publicado** (búsqueda sobre du y dv enteros y uniformes):

| Nivel (objetivo) | FIPS 203 | E8 | **Leech** |
|---|---|---|---|
| 512 (2⁻¹³⁹) | 768 B | 704 B (−8,3%; du=9, dv=4; 2⁻¹⁴²) | **696 B (−9,4%; du=9, dv=5; 2⁻¹⁶⁸)** |
| 768 (2⁻¹⁶⁴) | 1088 B | 1024 B (−5,9%; du=9, dv=5; 2⁻¹⁹⁴) | **984 B (−9,6%; du=9, dv=5; 2⁻¹⁸³)** |
| 1024 (2⁻¹⁷⁴) | 1568 B | 1408 B (−10,2%; du=10, dv=4; 2⁻²¹⁵) | **1376 B (−12,2%; du=10, dv=4; 2⁻¹⁸⁰)** |

**Validación de independencia en bloques de 24** (`leech_correlation_check.py`): las formas lineales dominantes (3w + Σ₂₃ ±w para Leech y Σ₈ ±w para el octad) se midieron con cifrados reales frente a coeficientes independientes. Cociente real/iid entre 0,97 y 1,13 con 1,2 a 2,4 M muestras, dentro del error estadístico. El supuesto se sostiene también para Leech.

**Comparación con Liu & Sakzad (2023, Leech24 + BW16):** usan otra métrica (más bits de mensaje con el mismo ciphertext, y la reducción la expresan en CER) y la aproximación gaussiana. En su configuración, ML-KEM-512/768 no alcanzan el objetivo de fallo (2⁻¹¹¹ y 2⁻¹³¹). Nuestro enfoque (clave fija de 256 bits, ciphertext más chico, cálculo exacto) cumple el objetivo en los tres niveles.

**Decodificador Leech de tiempo constante en C** (`v5/leech.c`, `test_leech.py`, `bench_leech.c`):
- **Construcción:** Golay vía Turyn, c = (a⊕x, b⊕x, a⊕b⊕x) con a, b ∈ H y x ∈ H′ (dos Hamming [8,4,4]). Verificado: 4.096 palabras con pesos (1, 759, 2576, 759, 1). Λ₂₄ = ∪_c (2c + 4D₂₄) ∪ (1 + 2c + 4(D₂₄+e₁)).
- **Decodificación de máxima verosimilitud exacta:** tablas por octeto (costo, paridad, flip más barato) y recorrido de los 2 × 4.096 candidatos con índices públicos. Selección, lookups y corrección de paridad por máscaras.
- **Corrección:** ida y vuelta OK (3.000 mensajes de 36 bits). Coincide con la búsqueda exhaustiva en numpy en 1.500 bloques (siempre encuentra el punto óptimo). Con ruido ML-KEM real (k=2, η=11, 600.000 bloques): fallo 1,27e-4 (76 eventos) contra 1,05e-4 predicho (cociente 1,21, unas 1,7σ; menos de 0,3 bit).
- **Velocidad:** unos **850.000 ciclos por mensaje (≈3.100× `poly_tomsg`)**. Impracticable tal como está. El decodificador óptimo de la literatura (Vardy–Be'ery, unas 3.600 operaciones por bloque, citado por Liu & Sakzad) bajaría a unas 29.000 operaciones por mensaje: sigue siendo más de 20× el decodificador E8.
- **Tiempo constante, hallazgo importante:** la primera versión **falló** dudect (|t| = 11,3). gcc -O3 había convertido las selecciones con máscara en saltos condicionales que dependen del secreto (`je`/`cmove` según la mitad y el bit ganadores). Se corrigió con una barrera de valor (`asm volatile`, la misma técnica que la referencia de Kyber). Después de la corrección: |t| = 2,16 y 2,60 (300.000 mediciones por corrida), y en el ensamblador quedan solo saltos de bucle.

**Pendiente o en contra de Leech:**
- Un decodificador de Leech en tiempo constante es mucho más complejo que el de E8. Liu & Sakzad reportan unas 3.595 operaciones en el peor caso y +34,5% de overhead, sin mediciones de side-channel.
- Hay que etiquetar los 2³⁶ puntos con bits. Sobran 32 bits de capacidad (288 contra 256).
- La búsqueda usó du y dv uniformes; con compresiones mixtas por componente el ranking podría cambiar.

## 9. v7 — Prototipo completo: ML-KEM (FIPS 203) contra ML-KEM+E8 (`v7/`)

**Implementación propia en C** del estándar completo: SHA3-256, SHA3-512, SHAKE128 y SHAKE256 (`fips202.c`); NTT, muestreo, CBD, compresión, K-PKE y ML-KEM KeyGen/Encaps/Decaps con transformación Fujisaki-Okamoto y rechazo implícito (`mlkem.c`). La variante E8 comparte todo el código salvo la codificación y decodificación del mensaje (decodificador AVX2 de `v4/e8fast.c`) y los valores de du y dv del modo compacto: 512 → (9,4), 768 → (9,5), 1024 → (10,4).

**Validación:**
- SHA-3/SHAKE: idénticos a `hashlib` de Python en tamaños de 0 a 2.000 bytes.
- **Vectores oficiales NIST ACVP (FIPS 203): 180/180 correctos.** keyGen 25/25, encapsulation 25/25 y decapsulation 10/10 (incluye ciphertexts inválidos) en los tres niveles.
- Variante E8: 20.000 rondas encapsular/desencapsular por nivel sin fallos, y rechazo implícito correcto ante ciphertexts alterados (200/200).

**Tamaños (bytes):**

| Nivel | Clave pública | Clave privada | Ciphertext estándar | **Ciphertext E8** | Diferencia |
|---|---|---|---|---|---|
| 512 | 800 | 1632 | 768 | **704** | −8,3% |
| 768 | 1184 | 2400 | 1088 | **1024** | −5,9% |
| 1024 | 1568 | 3168 | 1568 | **1408** | −10,2% |

Las claves no cambian. La seguridad no cambia (mismo ruido y mismo módulo). La tasa de fallo cumple el objetivo de NIST (sección 8: 2⁻¹⁴², 2⁻¹⁹⁴, 2⁻²¹⁵).

**Ciclos** (mediana de 5.000 repeticiones, estándar y E8 alternados en el mismo bucle):

| Nivel | KeyGen | Encaps | Decaps estándar | Decaps E8 | Diferencia en Decaps |
|---|---|---|---|---|---|
| 512 | igual (−0,1%) | −0,4% | 132.180 | 131.934 | −0,2% |
| 768 | igual (+0,1%) | 0,0% | 207.522 | 205.082 | −1,2% |
| 1024 | igual (0,0%) | −0,3% | 304.634 | 301.454 | −1,0% |

**Interpretación:** el decodificador E8 agrega unos 1.100 ciclos, pero el ciphertext más corto ahorra trabajo en la desencapsulación: el hash J(z‖c) absorbe un bloque Keccak menos en 768 y 1024, y la comparación de re-cifrado es más corta. **El efecto neto es nulo o levemente favorable.** Una primera medición con empaquetado bit a bit (lento) favorecía artificialmente a E8; se corrigió con empaquetado por palabras antes de medir.

**Límite:** esta implementación es de referencia, no optimizada como el código AVX2 oficial, así que los ciclos absolutos son más altos. En una implementación optimizada la desencapsulación es varias veces más rápida, y los ~1.100 ciclos del decodificador pesarían proporcionalmente más (del orden de pocos por ciento), compensados en parte por el ahorro de hashing. Queda por medir sobre el código optimizado oficial.

**Tiempo constante (desencapsulación completa, ciphertext válido fijo contra aleatorio, 50.000 mediciones):** |t| entre 0,96 y 1,78 en las seis variantes. No se detectan fugas.

**No implementado:** una API pública con aleatoriedad del sistema (los tests usan las funciones *internal* de FIPS 203 con semillas explícitas). Las verificaciones de entrada se agregaron en la sección 9.1.

### 9.1 Versión optimizada y verificaciones de entrada (`v7/mlkem_opt.c`)

**Optimización** con las mismas técnicas que el código de referencia oficial: Keccak-f[1600] completamente desenrollado (generado por script y verificado contra `hashlib`), absorción y extracción por palabras de 64 bits, NTT con multiplicación Montgomery y reducción Barrett, CBD por palabras, y para el mensaje estándar la fórmula exacta de `poly_tomsg` de la referencia.

| Bloque | Antes | Optimizado |
|---|---|---|
| Keccak-f (1 permutación) | 1.804 | 890 |
| SampleNTT (1 polinomio) | 16.474 | 6.382 |
| CBD η=2 | 8.238 | 1.414 |
| NTT | 6.342 | 1.764 |
| Desencapsulación ML-KEM-512 | 132.180 | 63.448 |

**Validación de la versión optimizada:**
- Vectores NIST ACVP: **180/180**.
- **Verificaciones de entrada FIPS 203** (módulo de ek, §7.2, y hash de dk, §7.3) implementadas en las dos versiones: **60/60 vectores ACVP keyCheck** correctos.
- **Equivalencia:** referencia y optimizada dan resultados idénticos byte a byte en las 6 variantes (1.200/1.200: ek, dk, c, K).

**Comparación final** (versión optimizada, mediana de 10.000 repeticiones alternadas; decodificación estándar con la fórmula oficial `poly_tomsg`):

| Nivel | KeyGen | Encaps | Decaps estándar | Decaps E8 | Diferencia en Decaps |
|---|---|---|---|---|---|
| 512 | +0,1% | −0,9% | 63.448 | 64.026 | **+0,9% (+578 ciclos)** |
| 768 | 0,0% | −0,4% | 98.116 | 98.052 | −0,1% |
| 1024 | −0,2% | −0,5% | 152.622 | 152.688 | 0,0% |

Tiempo constante en la desencapsulación completa (100.000 mediciones): |t| entre 0,83 y 3,03, sin fugas detectadas.

**Lección metodológica (dos veces):** el costo de E8 aparece o desaparece según lo eficiente que sea el código *estándar* contra el que se compara. Primero el empaquetado bit a bit, y después la decodificación estándar con división, favorecían artificialmente a E8. La comparación válida usa la fórmula oficial `poly_tomsg` (~270 ciclos) contra el decodificador E8 (~1.300 ciclos): **sobrecosto bruto de ~1.040 ciclos**, compensado en parte por hashear y comparar un ciphertext más corto.

**Estimación para el código oficial AVX2 (no medido acá):** el sobrecosto absoluto se mantiene en ~1.040 ciclos (el decodificador ya es AVX2), menos el ahorro de hashing. Para medirlo: compilar el benchmark oficial (por ejemplo, pq-crystals/kyber `avx2`, `test_speed`), anotar los ciclos de desencapsulación y sumar los ~1.040 ciclos. Ese porcentaje es la cota superior del costo de E8 en el código más rápido.

## 10. Modo de máxima seguridad: usar la ganancia de E8 para subir la seguridad (`v8/max_security.py`)

Objetivo original del framework: **superar la solidez de seguridad del estándar**. En lugar de achicar el ciphertext, la ganancia de E8 se usa para meter **más ruido** (η mayor) **manteniendo exactamente los tamaños de FIPS 203** (mismos k, du, dv) y la tasa de fallo exigida. Seguridad calculada con el estimador validado contra Kyber; fallo con el cálculo exacto.

| Nivel | Configuración | η1 | η2 | log₂ fallo | β | Bits clásicos | Bits cuánticos |
|---|---|---|---|---|---|---|---|
| 512 | FIPS 203 | 3 | 2 | −139,1 | 406 | 118,6 | 107,6 |
| 512 | **AFA-12 (E8)** | **4** | **3** | −146,6 | **424** | **123,8 (+5,2)** | **112,4 (+4,8)** |
| 768 | FIPS 203 | 2 | 2 | −165,3 | 624 | 182,2 | 165,4 |
| 768 | **AFA-12 (E8)** | **3** | 2 | −166,8 | **660** | **192,7 (+10,5)** | **174,9 (+9,5)** |
| 1024 | FIPS 203 | 2 | 2 | −175,2 | 874 | 255,2 | 231,6 |
| 1024 | **AFA-12 (E8)** | **3** | 2 | −189,4 | **923** | **269,5 (+14,3)** | **244,6 (+13,0)** |

- Mismos tamaños de claves y ciphertext, misma velocidad (como máximo +0,9% en desencapsulación). **Más margen ante mejoras futuras de los algoritmos de ataque, clásicos o cuánticos:** +10 bits equivalen a unas 1.000 veces más trabajo para el atacante.
- Es conservador: solo cuenta el ruido de la clave (η1). El η2 mayor en 512 también endurece el ciphertext.
- Límite honesto: el supuesto de dureza es el mismo (Module-LWE). Es más margen dentro de la misma familia, no un tipo nuevo de seguridad.
- Implementado y medido en la sección 10.3.
- **Corregido en la sección 14:** al estimar también el ataque dual y la instancia LWE del ciphertext, las ganancias quedan en +5,3 / +10,5 / +10,8 bits contando el ruido de redondeo, o en +6,1 / +5,3 / +7,3 sin contarlo. El +14,3 de 1024 era optimista.

### 10.1 Revisión de las 6 leyes con el objetivo "más seguridad"

| Ley (diagrama V6) | Diseño original | Resultado medido | Configuración final AFA-12 |
|---|---|---|---|
| 1. Materia y Fuerza (codificación) | 1 bit por coeficiente | **E8: el único punto que supera al estándar** | Codificación E8: +5/+10/+14 bits o fallo 2⁻²⁴¹ |
| 2. Canal Espacial (FO, r) | SHA3-512 + SHAKE-256 | Estándar y correcto. Una tasa de fallo menor ajusta la cota IND-CCA2 de la transformación FO | Igual a FIPS 203; r de CBD chico (no uniforme mod q) |
| 3. Masa de Higgs (ruido) | Gaussiana sobre E8 | **Pierde seguridad** (−17 bits a igual varianza; ataques reales lo confirman) | Ruido CBD con **η mayor**, posible gracias a la Ley 1 |
| 4. Restricción Lumínica (padding) | ISO 7816-4 | Correcto: oculta la longitud | Se mantiene |
| 5. Cristal E8 (matriz A) | A uniforme por SHAKE-128 | Correcto (validado con NIST) | Igual a FIPS 203 |
| 6. Fermión Fantasma (secreto) | Secreto sobre E8 | Pierde seguridad (mismo teorema que la Ley 3) | CBD con η mayor (igual que la Ley 3) |
| Capa DEM | K para SHAKE y para HMAC | **Riesgo:** reúso de clave | Derivar K_enc y K_mac separadas, o usar un AEAD estándar |

### 10.2 Robustez adicional para los próximos años (práctica de la industria, no medida acá)

- **Modo híbrido:** combinar el KEM con X25519, como hace TLS hoy con X25519MLKEM768. Si algún día se debilita la criptografía de lattices, sigue protegiendo la parte clásica, y al revés.
- **Agilidad criptográfica:** parámetros versionados para poder subir η o el nivel sin rediseñar.
- **Canales laterales de potencia y EM**, incluidos los ataques con redes neuronales: el decodificador E8 debe pasar esas pruebas en hardware real.

### 10.3 Modo de máxima seguridad y DEM corregida implementados (`v7/`)

**Modo de máxima seguridad** (`mlkem_opt.c`, ids 6 a 8: AFA12-512/768/1024-MaxSec):
- CBD con η=4 (extensión de AFA-12; FIPS 203 solo define η = 2 y 3), verificada contra la definición bit a bit de FIPS 203 junto con η=2 y 3 en 900 polinomios. Varianza empírica 2,000, teórica 2,0.
- Tamaños idénticos a FIPS 203. 20.000 rondas por nivel sin errores; rechazo implícito 200/200.
- Tiempo constante en la desencapsulación completa (100.000 mediciones): |t| = 2,06 / 1,75 / 1,19.
- Los vectores NIST de los modos estándar siguen pasando (180/180).

| Nivel | KeyGen | Encaps | Decaps | Ganancia de seguridad |
|---|---|---|---|---|
| 512 | +0,3% | +4,8% | +5,1% | +5,2 bits (clásico) / +4,8 (cuántico) |
| 768 | +6,4% | +2,7% | +2,9% | +10,5 / +9,5 |
| 1024 | +5,6% | +2,5% | +2,2% | +14,3 / +13,0 |

El costo es mayor que en el modo compacto porque más ruido requiere más salida de SHAKE (η=3 o 4 consume 2 bloques Keccak por polinomio en vez de 1), además del decodificador E8.

**Capa DEM corregida** (`dem.c`):
- **Separación de claves:** (K_enc‖K_mac) = SHAKE256("AFA12-DEM-KDF-v1"‖K). Se elimina el reúso de K del diseño V6.
- **Ley 4:** padding ISO/IEC 7816-4 a múltiplos de 64 bytes, que oculta la longitud.
- **Cifrado:** keystream SHAKE256 con separación de dominio y nonce de 16 bytes. **Autenticación:** HMAC-SHA3-256 sobre nonce‖C. **El tag se verifica en tiempo constante antes de descifrar** (encrypt-then-MAC, sin oráculo de padding).
- **Pruebas:** HMAC-SHA3-256 idéntico a `hmac` de Python (35 casos); ida y vuelta correcta con mensajes de 0 a 129, 500 y 1.000 bytes; **528/528 alteraciones rechazadas** (bit cambiado en C, tag, nonce o clave).

**Demo de punta a punta** (`demo_afa12.c`, flujo del diagrama V6 con AFA12-768-MaxSec): verificación de claves → encapsulación → DEM → desencapsulación → mensaje recuperado. Un ciphertext KEM alterado produce una clave no relacionada por el rechazo implícito, y la DEM rechaza. Un C alterado se rechaza antes de descifrar.

## 11. Cierre: configuración final de AFA-12

| Componente | Configuración final | Evidencia |
|---|---|---|
| Problema base | Module-LWE de ML-KEM (FIPS 203), q=3329, n=256, k=2/3/4 | 240/240 vectores NIST |
| Ley 1: codificación del mensaje | **E8** (8 bits cada 8 coeficientes, con la **distribución de bloques v13** de distancias distintas), decodificador AVX2 de tiempo constante | Fallo exacto validado con cifrados reales; dudect |
| Leyes 3 y 6: ruido y secreto | CBD estándar con **η mayor** gracias a la Ley 1 | Teorema y ataques reales: el ruido geométrico resta seguridad |
| Leyes 2 y 5: FO y matriz | Igual a FIPS 203 | Vectores NIST |
| Ley 4 + DEM | Padding ISO 7816-4, claves separadas, encrypt-then-MAC con HMAC-SHA3-256 | 528/528 alteraciones rechazadas |
| **Modo MaxSec** | Mismos tamaños, **+5 / +10 / +11 bits** contra el ataque más débil (o +6 / +5 / +7 sin contar el redondeo), +2% a +6% de tiempo | Secciones 10, 10.3 y 14 |
| Modo compacto | Ciphertext −8,3 / −5,9 / −10,2%, misma seguridad, ≤ +0,9% de tiempo | Sección 9.1 |

**Lo que AFA-12 no es:** no es más resistente que ML-KEM ante un ataque que rompa los lattices en general; el supuesto de dureza es el mismo. Para esa contingencia, la defensa de la industria es el modo híbrido con X25519.

**Pendientes para producción o publicación:** pruebas de canal lateral de potencia y EM (incluidos ataques con redes neuronales) sobre el decodificador E8; medición sobre el código oficial AVX2; modo híbrido con X25519; API con aleatoriedad del sistema; revisión externa por criptógrafos.

## 12. Camino a producción (`v9/`)

### 12.1 Modo híbrido con X25519
- **X25519 propio (RFC 7748)** con escalera de Montgomery de tiempo constante. Validación: vectores RFC 7748 (fijos y 1/1.000 iteraciones), **2.000/2.000 coincidencias con OpenSSL**. Al validar apareció un valor de vector mal transcrito de memoria; OpenSSL confirmó la implementación.
- **KEM híbrido estilo X-Wing:** K = SHA3-256("AFA12-HYBRID-v1" ‖ K_pq ‖ ss_x ‖ ct_x ‖ pk_x). La clave sigue protegida si cualquiera de los dos componentes resiste.
  - **Modo C:** ML-KEM-768 (FIPS 203) + X25519, compatible con certificación.
  - **Modo S:** AFA12-768-MaxSec + X25519, máxima seguridad.
- **API de producción:** aleatoriedad del sistema operativo (`getrandom`), verificaciones de entrada FIPS 203 obligatorias y borrado de secretos de memoria.
- **Pruebas** (ambos modos): combinador recalculado de forma independiente con OpenSSL y hashlib (300/300); alteración de la parte de lattices o de X25519 → clave distinta (100/100 cada una); claves inválidas rechazadas.

| Modo | pk | sk | ct | KeyGen | Encaps | Decaps | Tiempo constante |
|---|---|---|---|---|---|---|---|
| C: ML-KEM-768 + X25519 | 1.216 | 2.464 | 1.120 | 205.206 | 328.388 | 222.988 | \|t\| = 2,38 |
| S: AFA12-768-MaxSec + X25519 | 1.216 | 2.464 | 1.120 | 210.456 (+2,6%) | 330.636 (+0,7%) | 226.292 (+1,5%) | \|t\| = 1,97 |

(mediciones alternadas; X25519 ≈ 120.000 ciclos por operación y domina el costo)

### 12.2 Canales laterales de potencia (simulación, `power_sim.py`)
Metodología previa al silicio: fuga = peso de Hamming de cada valor intermedio + ruido gaussiano. Ataque perfilado (LDA, plantillas) que recupera bits del mensaje con **una sola traza**. Se compara la decodificación estándar (`poly_tomsg`) contra la E8 en condiciones idénticas.

| Ruido σ | Estándar: bits / bytes | **E8: bits / bytes** |
|---|---|---|
| 0,5 | 93,4% / 59,6% | 89,6% / 41,0% |
| 1 | 80,8% / 18,6% | 86,8% / 31,4% |
| 2 | 69,7% / 5,5% | 82,5% / 20,7% |
| 4 | 61,6% / 1,9% | **78,0% / 15,6%** |
| 8 | 56,1% / 1,1% | **74,9% / 15,2%** |

(azar: 50% / 0,39%)

**Hallazgo negativo importante:** el decodificador E8 hace más operaciones con el secreto y, en el modelo, **filtra bastante más por consumo eléctrico** que el estándar. Ninguno de los dos, sin enmascarar, resiste a un atacante con acceso físico.

**Regla de despliegue:**
- **Servidores, nube y equipos sin acceso físico del atacante:** el canal relevante es el tiempo, que ya está verificado. Modo S apto.
- **Tarjetas, IoT, HSM expuestos:** hace falta una implementación enmascarada de toda la desencapsulación. Es un proyecto aparte; mientras tanto, usar el Modo C.
- Falta confirmarlo en hardware real (por ejemplo, ChipWhisperer sobre Cortex-M4).

### 12.3 Auditoría externa: paquete listo (`v9/audit/`)
- **`AUDIT_PACKAGE.md`** (en inglés, para auditores internacionales): especificación exacta de las diferencias con FIPS 203, parámetros, híbrido y DEM; tabla de afirmaciones con su evidencia y fuerza; riesgos conocidos; comandos de reproducción; ruta de certificación; 7 preguntas concretas para el auditor.
- **`AFA12_KAT.json`:** 45 vectores KEM (las 9 variantes) + 10 híbridos, deterministas (SHA3-256 del archivo: 107eda61…8ed1). Las 15 entradas estándar coinciden con la versión validada por NIST.

### 12.4 Certificación FIPS
- AFA-12 **no es certificable** tal cual: no es un algoritmo aprobado.
- **Ruta compatible hoy:** Modo C (ML-KEM-768 + X25519), el mismo patrón que TLS (X25519MLKEM768). Para derivar claves que combinan un secreto de un esquema aprobado con otros secretos, NIST tiene la SP 800-56C Rev. 2, y la SP 800-227 trata combinadores de KEM. El formato exacto lo valida el laboratorio.
- **Cumplimiento + margen AFA-12:** híbrido triple ML-KEM-768 + AFA12-MaxSec + X25519, a un costo de unos 2,1 KB más por intercambio.
- **Largo plazo:** borrador IETF CFRG y comentario público para una revisión futura de FIPS 203, después de criptoanálisis independiente.

### 12.5 Versión enmascarada del decodificador E8 (preparada, `v10/`)

**Reformulación sin multiplicaciones** (`linear_check.py`): en E8 el error del candidato D8+½ es exactamente q menos el del candidato D8 (|e_B| = q − |e_A|), y la penalización por paridad es lineal en |e|. La comparación de distancias al cuadrado se reduce a
D = Σaᵢ − 4q + parA·(2q − 2·a_max) − parB·2·a_min, y se elige B si D > 0.
Solo usa sumas, comparaciones y selecciones. **Idéntica al decodificador validado en 400.000/400.000 bloques**, incluidos los empates.

**Decodificador enmascarado de primer orden** (`e8_masked.c`):
- **Entrada:** cuotas aritméticas w = w₀ + w₁ mod q (salida natural de una desencapsulación enmascarada). **Salida:** cuotas XOR del byte del mensaje (m = m₀ ⊕ m₁). El secreto nunca aparece sin repartir.
- **Gadgets:** AND seguro ISW (aleatoriedad nueva por operación); suma segura con acarreo por AND seguros (20 bits); conversión aritmética→booleana con resta condicional enmascarada de q; selecciones y extracción de signo por máscaras.
- **Opción de mezclado:** los 32 bloques se procesan en orden aleatorio.
- **Exactitud:** 96.000/96.000 bloques idénticos al decodificador validado, con y sin mezclado. Cada cuota de salida es uniforme (256 valores, media 127,3).

**Simulación de consumo** (`power_masked.py`; mismo modelo y ataque que la 12.2; σ = 1, el peor caso para el defensor):

| Prueba | Sin enmascarar | Enmascarado |
|---|---|---|
| TVLA primer orden (4.978 puntos, 20.000 + 20.000 trazas) | máx \|t\| = 86,0 | **máx \|t\| = 3,62: sin fuga** |
| Ataque LDA con 1 traza: bits / bytes | 86,8% / 31,6% | **50,1% / 0,42% (azar: 50% / 0,39%)** |
| TVLA segundo orden (producto de las 2 cuotas) | — | máx \|t\| = 78,9 (fuga esperable) |

**Interpretación:** el enmascarado de primer orden elimina los ataques de primer orden y de una sola traza en el modelo. Un atacante de segundo orden, que combina las dos cuotas, sigue teniendo señal, pero necesita muchas más trazas, y la cantidad crece con el ruido real. Para dispositivos de alta exposición se necesitan órdenes superiores (3 o más cuotas, mismos gadgets generalizados).

**Costo** (mediana, mismo banco de pruebas):

| Decodificador | Ciclos por mensaje |
|---|---|
| `poly_tomsg` ML-KEM (sin enmascarar) | 696 |
| E8 AVX2 (sin enmascarar) | 3.542 |
| **E8 enmascarado, primer orden** | **421.938** |
| E8 enmascarado + mezclado | 422.134 |

Es un prototipo de referencia, no optimizado. Camino conocido: **bitslicing** (32 coeficientes por palabra: cada AND seguro procesa 32 carriles a la vez, del orden de 30× menos) y sumador Kogge-Stone (profundidad logarítmica).

**Qué falta para una desencapsulación enmascarada completa:**
1. Enmascarar el resto del camino FO: NTT y productos (lineales, fáciles), compresión y CBD del re-cifrado (no lineales), Keccak enmascarado, comparación de ciphertexts enmascarada. Se puede partir de una implementación enmascarada publicada de ML-KEM y reemplazar solo el decodificador.
2. Generador de aleatoriedad de producción (el prototipo usa xorshift).
3. Validación en hardware real (ChipWhisperer / Cortex-M4): el modelo de valor por cuota no captura glitches ni transiciones entre registros, que en hardware pueden romper el enmascarado si no se separan las cuotas con cuidado.

### 12.6 Desencapsulación enmascarada completa y decodificador bitsliced (`v11/`)

**Punto 2: decodificador E8 enmascarado con bitslicing** (`masked_e8bs.c`, `masked_core.c`): mismo algoritmo sin multiplicaciones de la 12.5, pero cada palabra enmascarada lleva 32 carriles (un bloque de 8 coeficientes por carril). Así cada AND seguro procesa el mensaje completo de una vez.

| Decodificador | Ciclos por mensaje |
|---|---|
| E8 AVX2 sin enmascarar | 4.186 |
| E8 enmascarado v10 (bloque a bloque) | 421.922 |
| **E8 enmascarado v11 (bitsliced)** | **34.256 (12,3× más rápido)** |

Exactitud: 128.000/128.000 bloques idénticos al decodificador validado.

**Punto 1: desencapsulación FO enmascarada completa** (`masked_decaps.c`, `masked_keccak.c`). La clave secreta se reparte en dos cuotas aritméticas al cargarla (en el prototipo). A partir de ahí, todo el camino que depende del secreto opera sobre cuotas:

| Etapa | Técnica |
|---|---|
| u·s y NTT | lineal: se aplica a cada cuota por separado |
| decodificación E8 | A2B mod q bitsliced + decodificador 12.6 |
| G = SHA3-512(m′‖h) y J = SHAKE256(z‖c) | Keccak enmascarado (chi con AND seguro); verificado 175/175 contra hashlib |
| re-cifrado: CBD | ruido muestreado sobre bits enmascarados y convertido a cuotas aritméticas (B2A mod q) |
| codificación E8 de m′ | selección enmascarada de {0, 832, 1665, 2497} |
| comparación c′ = c | intervalos públicos de descompresión + reducción AND enmascarada |
| selección de K | sobre cuotas; solo se desenmascara la clave elegida |

- **Corrección de seguridad durante el desarrollo:** la primera versión desenmascaraba K′ y K̄ antes de elegir. Eso permitía un oráculo de comprobación de texto plano por canal lateral. Ahora se selecciona sobre cuotas.
- **Público por diseño (práctica estándar):** la bandera de validez del ciphertext.
- **Exactitud:** 720/720 sobre los 6 juegos AFA-12 (Compact y MaxSec). Hay 60 ciphertexts válidos y 60 alterados por juego; tanto la clave como el rechazo implícito coinciden con la versión sin enmascarar.

**Costo** (mediana; enmascarado y sin enmascarar intercalados; 1.000 repeticiones):

| Juego | Sin enmascarar | Enmascarado | Factor |
|---|---|---|---|
| AFA12-512-MaxSec | 62.532 | 397.704 | ×6,4 |
| AFA12-768-MaxSec | 101.176 | 540.214 | ×5,3 |
| AFA12-1024-MaxSec | 142.336 | 665.316 | ×4,7 |

Es el mismo orden que las implementaciones de ML-KEM enmascaradas de primer orden publicadas (de ×3 a ×7 en software).

**TVLA sobre la desencapsulación completa.** Condiciones: AFA12-768-MaxSec; modelo de peso de Hamming por palabra más ruido gaussiano σ = 1; 43.956 puntos por traza.

| Corrida | Trazas por clase | 1.er orden máx \|t\| | Puntos \|t\|>4,0 | Punto del máximo | 2.º orden máx \|t\| |
|---|---|---|---|---|---|
| fijo-vs-aleatorio | 10.000 | 4,50 | — | — | 99,6 |
| control aleatorio-vs-aleatorio | 10.000 | 4,31 | — | — | 3,7 |
| fijo-vs-aleatorio, semilla 2 | 30.000 | 4,78 | — | — | 171,5 |
| fijo-vs-aleatorio, semilla 3 | 30.000 | 5,04 | 3 | 13.928 | 171,1 |
| fijo-vs-aleatorio, semilla 4 | 30.000 | 4,50 | 5 | 38.688 | 170,2 |

**Interpretación de primer orden: sin fuga detectable.**
- Con 43.956 puntos, el umbral corregido por cantidad de pruebas (Bonferroni, α = 0,05) es \|t\| = 4,87. El 4,5 clásico es para un punto.
- Con 3× más trazas, una fuga real habría subido el máximo de ~4,5 a ~7,8. Sin embargo, se quedó entre 4,5 y 5,0.
- El único punto que superó 4,87 (5,04) no se repite en otra corrida independiente: el máximo cae en puntos distintos (13.928 contra 38.688). El criterio de doble corrida de ISO 17825 exige que el mismo punto falle dos veces, y no ocurre.
- La cantidad de puntos por encima de 4,0 (3 y 5) es compatible con el azar: lo esperado es 2,8 y P(≥5) = 0,15.

**Segundo orden:** fuga clara (\|t\| ≈ 170). Es lo esperado y lo que garantiza un enmascarado de primer orden: un atacante que combina las dos cuotas en el mismo instante sigue teniendo señal.

**Qué sigue faltando para un dispositivo real:**
1. CSPRNG de producción para las máscaras (el prototipo usa xorshift128+).
2. Clave secreta guardada ya enmascarada, con refresco de cuotas en cada uso (hoy se reparte al cargarla).
3. Validación en hardware (ChipWhisperer / Cortex-M4): glitches y transiciones entre registros no están en el modelo.
4. Orden superior (3 o más cuotas) para dispositivos de alta exposición, y sumador Kogge-Stone para bajar la latencia.

## 14. Re-estimación de seguridad: primal y dual, clave y ciphertext (`v12/`)

**Por qué:** la ganancia del modo MaxSec (sección 10) se había calculado con un solo ataque (primal) sobre una sola instancia (la clave). Hacía falta medir también el ataque dual y la instancia LWE del ciphertext. Esa instancia tiene secreto r ~ CBD(η₁) y error e₁ ~ CBD(η₂), más el redondeo de la compresión d_u.

**Herramienta** (`coresvp.py`): metodología core-SVP de la especificación de Kyber y del NIST. Costo clásico 0,292·β, cuántico 0,265·β. Ataque primal uSVP (estimación 2016) y ataque dual con ventaja ε = 4·exp(−2π²τ²) y amortización de listas de criba 2^0,2075β. Ambos escalan el secreto cuando su varianza difiere de la del error. *No se pudo usar el lattice-estimator oficial: la descarga desde GitHub está bloqueada en este entorno.*

**Validación contra Kyber (especificación de la ronda 3):**

| | Primal (nuestro / publicado) | Dual (nuestro / publicado) |
|---|---|---|
| ML-KEM-512 | 118,6 / 118 | 117,0 / 117 |
| ML-KEM-768 | 181,9 / 182 | 180,2 / 181 |
| ML-KEM-1024 | 254,9 / 256 | 252,6 / 253 |

**Resultado** (`coresvp_output.txt`, `weakest_output.txt`): seguridad del eslabón más débil, el mínimo entre {primal, dual} × {clave, ciphertext}.

| Nivel | Criterio | ML-KEM clásico | MaxSec clásico | Ganancia clásica | Ganancia cuántica |
|---|---|---|---|---|---|
| 512 | conservador (sin redondeo) | 113,5 | 119,5 | **+6,1** | +5,5 |
| 512 | con redondeo d_u | 117,0 | 122,3 | **+5,3** | +4,8 |
| 768 | conservador | 180,2 | 185,5 | **+5,3** | +4,8 |
| 768 | con redondeo | 180,2 | 190,7 | **+10,5** | +9,5 |
| 1024 | conservador | 252,6 | 259,9 | **+7,3** | +6,6 |
| 1024 | con redondeo | 252,6 | 263,4 | **+10,8** | +9,8 |

**Lectura:**
- **La ganancia se mantiene en todos los niveles y criterios** (de +5 a +11 bits clásicos), pero es menor que la declarada para 1024 (+14,3). En 768 y 1024, MaxSec deja η₂ = 2, así que el **ciphertext pasa a ser el eslabón más débil**: solo se endurece su secreto (η₁), no su error.
- El criterio "con redondeo" es el que usa el equipo de Kyber para su propia instancia del ciphertext en Kyber-512. El conservador no le da crédito al ruido de compresión.

**¿Se puede re-balancear?** (`rebalance.py`): busqué en η₁ ∈ {2..5}, η₂ ∈ {2..4}, d_u ∈ {9..11} y d_v ∈ {3..5}, con ciphertext ≤ FIPS y tasa de fallo con E8 ≤ la de FIPS. **Los juegos MaxSec actuales ya son los óptimos** en los tres niveles. Subir η₂ en 768 o 1024 rompe la tasa de fallo exigida. Ganar más requeriría aceptar más fallo o un ciphertext más grande.

### 14.1 Confirmación con el lattice-estimator oficial (`official_estimate.py`, commit f31b86f)

Corre todos los ataques del estimador oficial: uSVP, BDD, BDD híbrido, BDD MITM híbrido, dual, y dual híbrido con el modelo de costo MATZOV. El arora-gb y el BKW se excluyen por no ser competitivos. Las cifras son en log₂ de operaciones de bit. El ataque más fuerte resultó siempre el **dual híbrido**, que nuestro estimador no cubría. Eslabón más débil (mínimo entre clave y ciphertext):

| Nivel | Modelo | ML-KEM sin / con redondeo | MaxSec sin / con redondeo | Ganancia sin / con redondeo |
|---|---|---|---|---|
| 512 | core-SVP | 112,1 / 115,5 | 118,0 / 121,5 | **+5,9 / +6,0** |
| 512 | completo (MATZOV) | 136,0 / 139,7 | 141,6 / 144,9 | **+5,6 / +5,2** |
| 768 | core-SVP | 174,3 / 174,3 | 179,3 / 184,9 | **+5,0 / +10,6** |
| 768 | completo | 196,4 / 196,4 | 201,2 / 207,2 | **+4,8 / +10,8** |
| 1024 | core-SVP | 241,8 / 241,8 | 247,4 / 251,3 | **+5,6 / +9,5** |
| 1024 | completo | 262,3 / 262,3 | 267,6 / 271,7 | **+5,3 / +9,4** |

**Conclusión:**
- Los ataques más nuevos bajan las cifras absolutas de los dos esquemas por igual. **La ganancia de AFA-12 MaxSec se mantiene: +5 bits en todos los niveles y criterios, y hasta +11 si se cuenta el redondeo.** Coincide con nuestro estimador con una diferencia de 1 bit o menos (sección 14).
- **Dato relevante para 512:** con el modelo completo y el redondeo contado, AFA12-512-MaxSec llega a 2¹⁴⁴·⁹. Supera la referencia de la categoría 1 del NIST (2¹⁴³, el costo en compuertas de romper AES-128), que ML-KEM-512 no alcanza en este modelo (2¹³⁹·⁷). El NIST considera que ML-KEM-512 igual cumple, porque suma el costo de la memoria. Sin contar el redondeo, MaxSec queda en 2¹⁴¹·⁶.

**No cubierto:** la analogía con la criba cuántica es un modelo, no una medición. Tampoco se modelan los ataques que dependen de los fallos de descifrado (punto 2 del plan).

### 14.2 Otro compilador (clang) y Windows (`local/`)

- **Hallazgo:** compilado con **clang 18 -O3**, el decodificador E8 de referencia (`v4/e8codec.c`, `e8_decode`) tenía **fuga de tiempo**: dudect dio \|t\| ≈ 1.700, contra ≈ 1 con gcc. Clang convertía en saltos condicionales la lógica de máscaras que depende del secreto.
- **Corrección:** barreras de valor (`vb64`) en `ct_lt`, `ct_select`, `ct_abs` y en las máscaras de paridad. Después: \|t\| = 1,78 con clang y 0,90 con gcc. Los resultados no cambian: la suite completa de exactitud da lo mismo con los dos compiladores.
- **Impacto:** el KEM usa `e8_decode_avx2`, que nunca estuvo afectado. Con clang dan \|t\| ≤ 1,1 tanto `e8_decode_avx2` como `e8_decode_fast`; Leech da 1,03 y la desencapsulación completa de los 6 juegos queda en ≤ 1,85. Igual confirma que **el tiempo constante hay que verificarlo por compilador**.
- **Windows:** binarios con MinGW gcc y clang (`local/`). La prueba de equivalencia `xplat_kat.c` cubre los 9 juegos KEM, el híbrido C y S, la DEM y los decodificadores E8, con 25 a 40 casos deterministas cada uno. Da el **mismo digest global en Linux gcc, Linux clang, Windows gcc y Windows clang** (`f1e53458…` con la distribución v12; las de Windows probadas bajo Wine).
### 14.3 Pruebas en AMD Ryzen 5 5600X (Zen 3, Windows 11), binarios MinGW gcc 13 y clang 18

| Prueba | gcc | clang |
|---|---|---|
| Equivalencia (digest global `f1e53458…`, distribución v12) | idéntico a Linux | idéntico a Linux |
| Exactitud de los decodificadores (6,4 M bloques) | 0 diferencias | 0 diferencias |
| dudect de los decodificadores E8 (10 M mediciones) | ≤ 1,77 | ≤ 1,62 |
| Control positivo con fuga | 49.831 / 83.460 | 47.521 / 104.177 |
| dudect de la desencapsulación, 12 juegos (1 M / 100 k) | 11/12; **512-MaxSec 7,23** | 12/12, ≤ 3,49 |
| Costo MaxSec en la desencapsulación | +2,6% / +3,0% / +3,7% | +3,6% / +4,6% / +7,4% |

**Investigación del 7,23** (`02_recheck.bat`, `03_diag.bat`, `04_magnitud.bat`):
1. **Se repite** en el mismo binario gcc: 15,7, 11,4 y 10,4 con 1 M de mediciones. No aparece en Intel Xeon (nube), ni con el mismo `.exe` bajo Wine en Xeon (0,82), ni con clang en el Zen 3.
2. **No hay operaciones que dependan del secreto:** no hay divisiones en el código del KEM (revisado el desensamblado) ni saltos o tablas indexadas por secreto.
3. **Por etapas** (otro binario gcc del mismo código): la desencapsulación completa queda limpia (≤ 2,26). Solo hay 3,7 ciclos de diferencia en `pke_decrypt` según el ciphertext, que es público. El efecto depende de cómo quedó el código máquina en cada binario.
4. **Magnitud** (mismo armado que el binario que falló, 2 corridas de 1 M):

| Juego | Fijo vs aleatorio | Fijo vs aleatorio válido | Válido vs inválido (sin clase fija) |
|---|---|---|---|
| ML-KEM-512 (FIPS 203) | t 8,2 / 5,0; −10 / −6 ciclos | 0,9 / 1,1 | 1,7 / 1,5 |
| AFA12-512-MaxSec | t 14,4 / 11,8; −7 / −5 ciclos | 19,3 / 11,3; −9 / −5 ciclos | 1,4 / 1,9 |

**Conclusión:**
- El efecto es que **repetir exactamente la misma entrada hace la desencapsulación entre 5 y 10 ciclos más rápida sobre 85.000** (menos del 0,012%). Es típico de los predictores y buffers de la microarquitectura, que "aprenden" una ejecución repetida.
- **No distingue ciphertexts válidos de inválidos.** Es la propiedad que protege el rechazo implícito de FO, y la prueba C da t ≤ 1,9.
- **También afecta a ML-KEM-512 estándar** con el mismo compilador. No es propio de AFA-12.
- Un atacante solo aprendería si un ciphertext se repite, algo que ya conoce porque es quien lo envía.
- **Clasificación:** artefacto microarquitectónico dependiente de la entrada pública, de impacto bajo. Queda registrado como limitación conocida. **Recomendación:** distribuir binarios compilados con clang, que pasan todo en Zen 3, y repetir dudect en cada CPU objetivo.

## 15. Ataques por fallos de descifrado y corrección de la distribución de bloques E8 (`v13/`)

**Pregunta:** ¿puede un atacante elegir mensajes cifrados que fallen a propósito (failure boosting, D'Anvers et al. 2019) y aprovechar los fallos para obtener la clave? ¿El decodificador E8, que mira 8 coeficientes juntos, le da una ventaja nueva?

**Modelo.** El error al descifrar es δ_i = ⟨S, C_i⟩ + G_i:
- S = (e, −s) es secreto. C = (r, e₁ + error de compresión de u) y G = e₂ + error de compresión de v son conocidos por el atacante, que los elige probando mensajes m.
- **Criterio 1, norma:** buscar mensajes cifrados con ‖C‖² grande. Afecta igual al estándar y a E8.
- **Criterio 2, correlación** (solo E8): Cov(δ_a, δ_b) = σ²·A(b−a), donde A es la autocorrelación negacíclica de C. El estándar decodifica cada coeficiente por separado, así que la correlación no le importa. E8 falla a lo largo de 240 direcciones u, con varianza σ²·uᵀMu: si el atacante elige C con A(d) grande en las distancias d dentro de un bloque, sube la tasa de fallo de ese mensaje cifrado.
- **Costo del ataque** para obtener el primer fallo: trabajo W = 1/(α·β) y consultas Q = 1/β, con el límite del NIST Q ≤ 2⁶⁴. α es la fracción de mensajes "débiles" y β su tasa de fallo.

**Herramientas:**
- `failure_boost.py`: criterio de norma, con cota de Chernoff para α.
- `corr_boost.py` / `corr_layout.py`: Monte Carlo de la correlación.
- `corr_chernoff.py`: cota de Chernoff para la forma cuadrática Z = CᵀQC/CᵀC con los autovalores negacíclicos μ_j = Σ_l w_l cos(π(2j+1)l/256). La cola resultó más pesada que una gaussiana, y esta cota queda unas 10 veces por encima del Monte Carlo, así que favorece al atacante.
- `attack_cost.py`: costo total del ataque.

**Validación con cifrados reales** (`validate_corr.py`: k=2, η=14 inflado, 80.000 mensajes cifrados): **299 fallos observados contra 303,1 predichos** por el modelo, mensaje por mensaje. Por grupos de correlación: 162/147, 103/123, 30/29 y 4/3,6.

**Resultado 1, criterio de norma solamente:** con Q ≤ 2⁶⁴ el trabajo es 2⁵¹⁵ o más en todos los juegos, incluido ML-KEM. Sin límite de consultas cuesta ≈ 1/DFR, igual para todos.

**Resultado 2, correlación: punto débil del diseño anterior.** Con bloques de **8 coeficientes consecutivos**, las distancias cortas se repiten dentro del bloque: la distancia 1 aparece 7 veces, la 2 seis veces, etc. Basta forzar una sola distancia para afectar muchos pares a la vez. El desvío de la dirección más desfavorable llega a 2,96 veces el de una sola correlación.

**Corrección (implementada):** una **distribución pública de bloques** (`v7/e8_layout.h`, buscada con `layout_search.py`), en la que cada bloque toma 8 posiciones cuyas 28 distancias internas, como clases min(d, 256−d), son todas distintas. Forzar una distancia ahora afecta un solo par por bloque, con peso ¼; el desvío máximo baja a 1,32 veces, y cada distancia se repite en a lo sumo 8 bloques. La distancia 128 tiene correlación cero por la estructura del anillo (A(256−d) = −A(d)).

Tasa de fallo que puede alcanzar el atacante después de invertir 2⁶⁴ o 2¹²⁸ de trabajo en elegir mensajes (log₂; cuanto más negativa, mejor):

| Juego | 2⁶⁴ consecutiva | 2⁶⁴ **nueva** | 2⁶⁴ ML-KEM | 2¹²⁸ consecutiva | 2¹²⁸ **nueva** | 2¹²⁸ ML-KEM |
|---|---|---|---|---|---|---|
| 512-Compact | −65,9 | **−81,4** | −103,7 | −56,6 | **−65,2** | −90,8 |
| 512-MaxSec | −69,9 | **−90,3** | −103,7 | −58,6 | **−71,0** | −90,8 |
| 768-Compact | −90,7 | **−116,5** | −127,9 | −77,2 | **−91,2** | −114,3 |
| 768-MaxSec | −89,4 | **−114,1** | −127,9 | −73,9 | **−93,3** | −114,3 |
| 1024-Compact | −121,4 | **−150,3** | −141,3 | −99,2 | **−125,0** | −127,2 |
| 1024-MaxSec | −106,6 | **−134,1** | −141,3 | −85,8 | **−113,4** | −127,2 |

**Costo total para obtener el primer fallo con Q ≤ 2⁶⁴** (log₂ de trabajo, cota conservadora):

| Juego | Solo norma | Consecutiva + correlación | **Nueva + correlación** | Seguridad del nivel (estimador oficial) |
|---|---|---|---|---|
| 512-Compact | 615 | **122 ⚠️** | **185** | ≈ 140 |
| 512-MaxSec | 582 | **143 ⚠️** | **225** | 144,9 |
| 768-Compact | 1544 | 244 | **358** | ≈ 196 |
| 768-MaxSec | 1269 | 237 | **379** | 207,2 |
| 1024-Compact | — | 549 | **783** | ≈ 262 |
| 1024-MaxSec | 1933 | 331 | **545** | 271,7 |

**Conclusión:**
- **Con la distribución consecutiva, el nivel 512 quedaba en el límite o por debajo de su seguridad declarada** frente a este ataque, siempre según la cota conservadora, que favorece al atacante.
- **Con la nueva distribución, todos los juegos superan su nivel con margen amplio:** 2¹⁸⁵ o más contra 2¹⁴³ en el 512.
- Es un punto débil propio de cualquier código de varias dimensiones (E8, Leech), no descripto en el trabajo previo (Liu & Sakzad). Hay que tenerlo en cuenta en cualquier propuesta de este tipo.

**Verificación de la corrección:**
- Tasa de fallo promedio sin cambios con cifrados reales (`dfr_layout_check.py`): 2,11·10⁻⁵ en las dos distribuciones (81/81 en 3,84 M de bloques). En otra semilla dieron 28 contra 32, con fallos en mensajes distintos (0 en común), lo que confirma que no era un error del test.
- NIST ACVP 180/180, controles de claves 60/60, referencia = optimizada (1.200/1.200).
- Desencapsulación enmascarada = normal (720/720). Híbrido OK.
- **Vectores de prueba de AFA-12 regenerados** (`v9/audit/AFA12_KAT.json`, SHA3-256 `c1ccef4b…05e0`): los 15 de ML-KEM y los 5 del híbrido C no cambian; cambian los 30 de E8 y los 5 del híbrido S. Los anteriores quedan en `AFA12_KAT_v12_consecutive.json`.
- Nuevo digest multiplataforma `d1174378…c3a5`, idéntico en Linux gcc y clang, y en Windows bajo Wine.

**Pendiente:**
- Análisis del boosting direccional después del primer fallo (D'Anvers et al. 2020), que abarata los fallos siguientes. Afecta a todos los esquemas; el primer fallo es el paso dominante.
- Cola exacta de la correlación más allá de la cota de Chernoff (por ejemplo, con muestreo por importancia).

## 16. Cierre de los análisis pendientes de fallos (`v14/`)

### 16.1 Estimación ajustada del ataque por correlación (`corr_saddle.py`)
La cota de Chernoff de la sección 15 se reemplazó por la **aproximación de punto silla de Lugannani-Rice** para la forma cuadrática.
- **Validación contra Monte Carlo** con las entradas reales (CBD más redondeo), en AFA12-512-MaxSec:

| Probabilidad simulada | 10⁻² | 10⁻³ | 10⁻⁴ |
|---|---|---|---|
| Dirección tipo A, distancia 1: punto silla | 1,02·10⁻² | 0,97·10⁻³ | 0,97·10⁻⁴ |
| Dirección tipo B, bloque v13: punto silla | 1,01·10⁻² | 1,04·10⁻³ | 1,15·10⁻⁴ |
| Chernoff (v13), para comparar | 6,8–7,2·10⁻² | 0,8–0,9·10⁻² | 1,0–1,3·10⁻³ |

- **Costo del primer fallo** con Q ≤ 2⁶⁴ (log₂), punto silla contra la cota de Chernoff de v13:

| Juego | Consecutiva | **Distribución v13** | v13 según Chernoff | Seguridad del nivel |
|---|---|---|---|---|
| 512-Compact | 127,4 | **192,0** | 185,1 | ≈ 140 |
| 512-MaxSec | 147,6 | **232,4** | 225,3 | 144,9 |
| 768-Compact | 249,0 | **365,4** | 358,3 | ≈ 196 |
| 768-MaxSec | 242,2 | **386,6** | 379,1 | 207,2 |
| 1024-Compact | 554,4 | **790,7** | 783,1 | ≈ 262 |
| 1024-MaxSec | 336,4 | **552,3** | 544,7 | 271,7 |

- **Lectura:** la estimación ajustada sube los costos unos 7 bits. **Las conclusiones de la sección 15 se mantienen:** con bloques consecutivos, el nivel 512-Compact quedaba por debajo de su seguridad (2¹²⁷); con la distribución v13, todos los niveles quedan muy por encima.

### 16.2 Ataque después del primer fallo: boosting direccional (`directional.py`)
Se usa el modelo geométrico de D'Anvers, Rossi y Virdia (EUROCRYPT 2020).
- Cada fallo obtenido revela un vector con coseno t respecto del secreto. Con n fallos, el atacante tiene una estimación con coseno ρₙ y busca mensajes cifrados alineados con ella.
- Supuestos a favor del atacante: conoce el coeficiente o la dirección que falló, y cada fallo aporta información completa.
- E8 tiene más direcciones posibles (32×240 contra 256), lo que le da al atacante más opciones de alinear. Pero su umbral es √2 veces más alto en proporción.

Trabajo (log₂) para obtener el fallo n+1, dados n fallos:

| Juego | t | n=1 | n=2 | n=3 | n=5 |
|---|---|---|---|---|---|
| ML-KEM-512 | 0,51 | 84,7 | 81,3 | 80,4 | 80,2 |
| AFA12-512-Compact | 0,54 | 145,4 | 126,4 | 121,1 | 117,1 |
| AFA12-512-MaxSec | 0,55 | 147,9 | 129,7 | 123,5 | 119,4 |
| ML-KEM-768 | 0,45 | 103,0 | 91,7 | 89,1 | 86,7 |
| AFA12-768-Compact | 0,46 | 302,2 | 231,6 | 209,1 | 193,1 |
| AFA12-768-MaxSec | 0,48 | 200,2 | 161,1 | 147,9 | 137,1 |
| ML-KEM-1024 | 0,37 | 280,6 | 207,7 | 184,4 | 164,7 |
| AFA12-1024-Compact | 0,48 | 284,7 | 214,8 | 192,0 | 174,5 |
| AFA12-1024-MaxSec | 0,39 | 360,4 | 258,8 | 225,7 | 201,1 |

**Lectura:**
1. Igual que en la literatura sobre Kyber, **una vez obtenido el primer fallo, los siguientes son mucho más baratos**, en todos los esquemas. Por eso la seguridad frente a este ataque la define el **primer fallo**.
2. **Los fallos siguientes son más caros en AFA-12 que en ML-KEM en todos los niveles:** 2¹⁴⁵ contra 2⁸⁵ en 512, y 2²⁰⁰ contra 2¹⁰³ en 768-MaxSec. El margen extra de E8 juega a favor en esta etapa.
3. Ataque completo = primer fallo (2¹⁹² o más con la distribución v13) más los siguientes. **Todos los juegos quedan por encima de su nivel.**
4. **Multi-objetivo:** si el atacante ataca 2ᴷ claves a la vez, el primer fallo se abarata hasta 2ᴷ. Es igual para ML-KEM, y los niveles del NIST se definen sin contarlo. Con la distribución v13, 512-Compact aguanta K ≤ 49 sin bajar de 2¹⁴³.

### 16.3 Código Leech: exposición a la correlación (`leech_layout.py`)
Medida de exposición: cuántas veces amplifica la peor dirección mínima la correlación de una sola distancia.

| Código | Bloques consecutivos | Mejor distribución encontrada |
|---|---|---|
| E8 (8 coordenadas) | 2,96 | **1,32** (v13; todas las distancias distintas) |
| Leech (24 coordenadas) | 4,90 | 2,14 (búsqueda local; con 276 pares y solo 127 clases de distancia, es imposible que sean todas distintas) |

- **Conclusión:** los códigos de más dimensiones están **más expuestos**. La defensa por distribución ayuda (4,90 → 2,14), pero no llega al nivel de E8.
- Es otro motivo, además de la velocidad (850 k ciclos), para preferir E8. Lo dejo como recomendación general: **cualquier código multidimensional usado para codificar el mensaje en LWE tiene que analizarse frente a este ataque.**

## 13. Próximos pasos

1. Integrar el codec en una implementación completa de ML-KEM y medir el costo relativo sobre la desencapsulación.
2. Repetir las pruebas de tiempo constante en otra CPU y otro compilador (clang, ARM).
3. Escribir la nota técnica: cifras exactas para ML-KEM-512/768/1024, reconciliación con Liu & Sakzad, validación con cifrados reales y el codec de tiempo constante con sus pruebas.
