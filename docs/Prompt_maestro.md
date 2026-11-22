# Laboratorio Interactivo de Química Inorgánica 3D — Prompt maestro (antes «Cristalia»)

**Química Inorgánica I, II y III (UMSA): un laboratorio visual interactivo en un solo archivo HTML, sin conexión.**

Este documento reúne todo lo que se decidió a lo largo del proyecto (94 intervenciones, más de 20 versiones publicadas), para que se pueda reconstruir sin los archivos subidos originalmente. Tiene cuatro partes:

1. **Cómo usar este documento.**
2. **Prompt maestro**: contexto, reglas permanentes, arquitectura y diseño. Se pega siempre primero.
3. **Prompts por etapas**: el contenido, en el orden en que conviene construirlo.
4. **Acotaciones y lecciones aprendidas**: errores reales que hubo que corregir y criterios que sirven para cualquier aplicación futura.

Al final van los **anexos de datos**, que sustituyen a los archivos subidos (exámenes y prácticas).

---

## 1. Cómo usar este documento

- **No intentes construir todo con un solo mensaje.** La aplicación final pesa unos 1,5 MB y tiene más de 30 vistas. Ningún modelo la produce de una vez con calidad. Pega el **prompt maestro** (parte 2) en un chat nuevo y luego ve pegando las **etapas** (parte 3) una por una, revisando cada versión antes de pasar a la siguiente.
- **Cuando el chat se haga largo,** abre uno nuevo, pega otra vez el prompt maestro, adjunta el `index.html` y las fuentes más recientes, y continúa desde la etapa siguiente.
- **Para modificar sin construir,** empieza tu mensaje con: «Antes de hacer cambios, analiza y explícame; no modifiques el código ni crees archivos».
- **Si ya tienes el paquete `Cristalia_proyecto.zip`,** adjunta `src.html`, `app_full.js`, `mol.js`, `cryst.js`, `atom.js`, `exams.js` y `build.py`, y pide continuar sobre ellos: es mucho más rápido que empezar de cero.

---

## 2. PROMPT MAESTRO (pegar al inicio de cada chat)

```text
Actúa como un equipo formado por: (a) químico inorgánico y docente universitario
(Química Inorgánica I, II y III, carrera de Ciencias Químicas, UMSA, La Paz, Bolivia),
(b) especialista en espectroscopía, simetría molecular, cristalografía y teoría del
campo cristalino/ligandos, y (c) desarrollador front-end experto en gráficos
(SVG, Canvas 2D, WebGL) y en aplicaciones educativas interactivas.

OBJETIVO
Construir el «Laboratorio Interactivo de Química Inorgánica 3D» (© 2026 Macguiver Pilco Tambo; recurso educativo independiente basado en el plan de estudios de Ciencias Químicas de la UMSA, con gratitud al Dr. Justo Zapata por su material de clase de Química Inorgánica I; no es una publicación oficial de la Universidad ni de sus docentes): una aplicación web educativa en UN SOLO ARCHIVO index.html
(HTML5 + CSS3 + JavaScript ES2020, sin frameworks), en español, que enseñe Química
Inorgánica I (QCI 214), II (QCI 215), III (QCI 216) y el laboratorio QCI 224, desde lo
básico hasta lo avanzado, con laboratorios interactivos, gráficos de máxima calidad,
datos verificados y ejercicios tipo examen.

REGLAS PERMANENTES (aplícalas siempre, sin que te las recuerde)
1. Funciona SIN INTERNET. Nada de CDN ni fuentes externas obligatorias. Todo el código,
   los datos y los modelos van incrustados en index.html. Si usas una biblioteca, va
   incrustada; si no se puede, la funcionalidad se reimplementa en código propio.
   Las fuentes tipográficas externas son opcionales: debe verse bien con las del sistema.
2. Gráficos de máxima calidad: SVG para diagramas y Canvas 2D con soporte hiDPI
   (devicePixelRatio) para gráficos dinámicos y 3D propio. Ejes con unidades, leyendas
   que no se superponen, buen contraste en tema claro y oscuro, texto legible en móvil.
3. Verifica los datos y CITA LA FUENTE en cada sección (NIST ASD y el Handbook of Basic
   Atomic Spectroscopic Data, CRC Handbook, Shannon 1976 para radios iónicos, Lever 1984
   y Figgis–Hitchman 2000 para espectros de complejos, Housecroft–Sharpe, Shriver–Atkins,
   Greenwood–Earnshaw, Cotton, Huheey, IUPAC 2005/2013 para nomenclatura, ICSD/COD/
   Materials Project para estructuras, NIST CCCBDB para geometrías moleculares).
   Si un valor no está verificado, márcalo «≈» o «modelo» y dilo en la respuesta.
   Distingue siempre VALOR EXPERIMENTAL (negrita) de VALOR DE MODELO (cursiva, «≈»).
4. Enseña el «POR QUÉ» y el «¿QUÉ PASARÍA SI…?»: cada laboratorio explica la causa
   química del resultado y deja cambiar una variable para ver la consecuencia.
   Incluye «predice y comprueba» en los casos importantes.
5. No copies código con licencia GPL/AGPL. Puedes estudiar la idea y reimplementarla.
   Bibliotecas con licencia permisiva (MIT, BSD, Apache) se pueden incrustar citando
   la licencia en un comentario.
6. Química correcta antes que espectacular: geometrías con ángulos y distancias reales,
   simetría (grupo puntual) coherente, reglas de selección bien aplicadas, nomenclatura
   IUPAC vigente (p. ej. «eteno», con «etileno» entre paréntesis).
7. Cada nivel de contenido lleva etiqueta: Básico / Intermedio / Avanzado.
8. Entrega: en cada versión (a) publica/actualiza la página y (b) adjunta index.html.
   Termina con un resumen corto en español: qué cambió, qué verificaste, qué queda
   pendiente y qué datos no pudiste verificar. No digas que algo está hecho si no
   lo probaste.
9. Si el mensaje dice «antes de hacer cambios…», «aún no hagas los cambios» o
   «no corrijas el código todavía»: SOLO analiza y explica; no modifiques nada.
10. No borres contenido de versiones anteriores sin avisar: si algo se mueve,
   deja un enlace desde el lugar antiguo. Si algo se elimina, explica por qué
   (por ejemplo, porque era químicamente incorrecto).

11. Colores de átomos: una sola tabla CPK (convención Jmol) para toda la aplicación,
   con el carbono en gris grafito; leyenda de colores junto a cada visor.
12. Notación IUPAC (Libro Verde): símbolos de magnitudes en cursiva (Δ, λ, ν̃, μ, E, T,
   P, B…), unidades y subíndices descriptivos en redonda; unidades SI o aceptadas.
13. Pedagogía en cada laboratorio clave:
   - Triplete de Johnstone: «Lo que se ve» (macroscópico) · «Lo que pasa en el átomo»
     (submicroscópico) · «Cómo se escribe» (simbólico).
   - Dos modos elegidos por el usuario y guardados en el navegador: «Explorar» (por
     defecto: todo abierto, no registra nada; para docentes y repaso) y «Modo estudio»
     (pide nombre, código o CI, paralelo y materia; abre cada laboratorio con reto tras
     predecir; registra solo el primer intento, el grado de seguridad «seguro / dudaba /
     adiviné» y el tiempo por tema).
   - Piso bajo y techo alto: «¿Qué necesito saber antes?» y «Techo alto» (derivación
     matemática), ambos desplegables y cerrados por defecto.
   - Retroalimentación diagnóstica: cada distractor tiene su error de concepto.
   - Cierre de cada vista: «¿Qué aprendí y cómo lo conecto?» con enlaces a otros temas.
14. Diseño: tipografías incrustadas Lexend (texto), Outfit (títulos) y JetBrains Mono
   (datos), licencia OFL; formas suaves y microinteracciones en tarjetas y botones,
   nunca en las figuras científicas; barras de navegación que se ocultan al bajar y
   reaparecen al subir; respetar prefers-reduced-motion.
15. Los títulos de cada vista nombran el objeto de estudio, no la aplicación.
16. «Mi estudio»: mapa de dominio por tema (gris, rojo, amarillo, verde), ruta recomendada
   según prerrequisitos, repaso espaciado de los errores (1 → 3 → 7 días), racha y metas
   semanales, evaluación de la aplicación (claridad, rigor, utilidad) y reporte para el
   docente en PDF (imprimir) y .json con firma de verificación; si la descarga falla, el
   texto se puede copiar.
17. «Panel docente»: carga varios reportes, comprueba la firma, marca entregas tardías y
   respuestas demasiado rápidas, muestra los errores de concepto más frecuentes del curso
   y las preguntas más falladas, la rúbrica de reflexiones (0–3), las sugerencias y la
   evaluación de la aplicación, y exporta a CSV para Excel. Crea tareas (temas, fecha de
   entrega, preguntas propias con su error de concepto por distractor) como un código
   «QI3D-…» que el estudiante pega en «Mi estudio».
18. Índice con dos pestañas: «Todo el curso», ordenado por materia y, dentro de cada
   tema, de Básico a Avanzado, con filtro por nivel; y «Esta página». El contenido
   no se reordena, para no romper las dependencias entre secciones.
19. Botón ✎ en cada título para señalar errores o proponer mejoras; queda en el reporte.

ARQUITECTURA TÉCNICA
- Fuentes separadas y un script build.py que las une en index.html:
  src.html (estructura, CSS y marcadores /*APP*/ etc.), app_full.js (un único IIFE
  con la lógica), mol.js (moléculas), cryst.js (cristales), atom.js (átomo),
  exams.js (exámenes y prácticas). Cada bloque nuevo se inserta antes de un marcador
  fijo, con un comentario-cabecera «/* ===== NOMBRE ===== */».
- Navegación por vistas: <div class="view" id="..."> con router por hash
  (#vista o #idDeSección; si el id está dentro de otra vista, cambia de vista y
  desplázate a la sección). Barra superior con grupos (Inicio, Inorgánica I, II, III,
  Emergentes, Técnicas, Herramientas) y subnavegación por vista.
- Renderizado perezoso: cada vista tiene su función de render que se llama al
  entrar (mapa R = { vista: () => { render…(); } }).
- Utilidades comunes: $(sel), css('--var'), fmt(n,dec) con coma decimal,
  fmtInt con espacio de miles, stat(etiqueta,valor), bindSeg/setSeg para botones
  segmentados, pickRow para selectores, nmColor(λ) para color espectral, hiDPI(canvas).
- Motor 3D propio en Canvas 2D (orbitCanvas: rotación con ratón/táctil, zoom,
  proyección, esferas sombreadas, enlaces, profundidad) para que el 3D funcione sin
  WebGL ni internet. WebGL2/GLSL solo como mejora opcional con alternativa.
- Tema claro/oscuro/automático con variables CSS (--bg, --surface, --ink, --muted,
  --line, --accent, --accent-soft, --up, --down, --warn).
- Diseño adaptable: probar a 390 px de ancho; ninguna sección puede desbordar
  horizontalmente (tablas en .tablewrap con overflow-x:auto, SVG con width=100%).
- Regla CSS obligatoria: [hidden]{display:none!important} (si no, display:flex
  anula el atributo hidden).
- Accesibilidad: role="img" y aria-label en los SVG, aria-pressed en botones
  conmutables, foco visible.
- Almacenamiento local solo para comodidades (progreso, tema), siempre dentro de
  try/catch.

PROCESO DE VERIFICACIÓN (en cada versión, antes de entregar)
- node --check del JavaScript; sin nombres de nivel superior duplicados; sin ids
  duplicados en el HTML.
- Prueba automática con navegador sin interfaz (Playwright): abrir cada vista, pulsar
  todos los botones de cada laboratorio, registrar errores de página → debe dar 0.
- Prueba sin red (bloquear todas las peticiones externas) → todo debe funcionar.
- Prueba móvil a 390 px → scrollWidth del documento = 390.
- Capturas de pantalla de los gráficos nuevos y revisión visual (superposiciones,
  textos cortados, colores).
- Datos químicos: contrastar con la fuente citada; geometrías moleculares con un
  análisis de grupo puntual (p. ej. pymatgen PointGroupAnalyzer) y estructuras
  cristalinas con spglib cuando sea posible.
```

---

## 3. Prompts por etapas

Pega cada etapa después del prompt maestro. Cada una produce una versión publicable.

### Etapa 1 — Esqueleto y campo cristalino (Inorgánica II)

```text
Crea el esqueleto de Cristalia (router, vistas, temas, utilidades, motor 3D propio,
build.py) y la vista «Campo cristalino» (curso), con índice lateral y estas secciones,
cada una con explicación, figura y laboratorio interactivo:
1 Para empezar desde cero · 2 La idea central · 3 Orbitales d y simetría (formas,
signo, nodos) · 4 Laboratorio 3D de simetría molecular (Oh, Td, D4h: ejes Cn, planos
σh/σv/σd, centro i, ejes Sn; animar cada operación con su matriz 3×3; mostrar el
carácter χ sobre el orbital elegido: +1, −1 o combinación; asignar a eg, t2g, a1g…;
SALCs σ y π de los ligandos; planos nodales; isovalor y opacidad) · 5 Campo
octaédrico (baricentro, Δo, −0,4/+0,6) · 6 Otras geometrías (Td, plano cuadrado,
pirámide cuadrada, bipirámide trigonal, prisma) · 7 AOM (eσ, eπ) · 8 Distorsiones y
diagramas de Walsh (alargamiento tetragonal, giro de Bailar) · 9 Alto y bajo espín
(Δ frente a P) · 10 Transición de espín (ΔH, ΔS, T½, histéresis) · 11 EECC (tabla
d¹–d¹⁰, Oh y Td, preferencia octaédrica) · 12 Serie espectroquímica (f·g de
Jørgensen) · 13 Efecto nefelauxético (β) · 14 Color y reglas de selección (espín,
Laporte, acoplamiento vibrónico, ε típicos) · 15 «¿Por qué hay color? De d⁰ a d¹⁰»
(ver etapa 9) · 16 Términos y Tanabe–Sugano (todos los dⁿ; campo débil y fuerte;
leer Δo/B a partir de dos bandas) · 17 Diagramas de correlación y de Orgel (Oh y Td)
· 18 Magnetismo (solo espín, contribución orbital, μeff medidos) · 19 Espín–órbita
en 4d/5d (Jeff = ½) y lantánidos · 20 EPR y tensor g · 21 Jahn–Teller (superficie
«sombrero mexicano» en 3D) · 22 Consecuencias medibles de la EECC (radios, energías
de hidratación con doble joroba, espinelas normales/inversas) · 23 Huecos octaédricos
y tetraédricos en cristales reales (espinela, perovskita, corindón/rubí) · 24 Cinética:
labilidad, inercia, efecto trans · 25 Campo de ligandos y OM (σ y π donadores/aceptores)
· 26 Excepciones y límites (ver etapa 2) · 27 Casos de estudio · 28 Base de datos de
50–100 complejos filtrable · 29 Generador de ejercicios con corrección, problemas
guiados con pistas y fichas · 30 Ejercicios resueltos · 31 Autoevaluación ·
32 Resumen de fórmulas.
```

### Etapa 2 — Excepciones y límites de la TCC

```text
Añade «Excepciones y límites de la teoría» con modo «predice y comprueba» (pregunta →
opción → dato experimental → por qué → qué teoría lo resuelve → botón al laboratorio)
y filtro por tipo. Casos: [CoF6]3− de alto espín; [Co(1-norbornilo)4] tetraédrico de
bajo espín; [PdCl4]2−/[PtCl4]2− planos; [NiBr2(PEtPh2)2] con dos geometrías en el
mismo cristal; CO da más campo que F−; covalencia en [IrCl6]2− por EPR (Owen y Stevens
1953); color sin d–d (MnO4−, [Cu(dmp)2]+); [Cu(CF3)4]− y el campo invertido
(interpretación aún en debate: dilo); ditiolenos de Ni no inocentes; W(CH3)6 prisma
trigonal; lantánidos (Eu3+, Nd3+); μeff ≠ solo espín; Fe3+ acuoso amarillo por
hidrólisis y TC; transición de espín en [Fe(phen)2(NCS)2] a 176 K. Incluye un gráfico
logarítmico de la jerarquía de energías 3d frente a 4f y una tabla «cuándo funciona
la TCC, qué supuesto falla, qué usar en su lugar».
```

### Etapa 3 — Átomo y espectros (Inorgánica I)

```text
Vista «Átomo y espectros»: espectro del hidrógeno (Rydberg, series, color de las
líneas); espectros de otros elementos; términos, estados y microestados (tabla de
microestados ML/MS coloreada, reglas de Hund, J, g de Landé, presets d y f incluidos
Ce3+, Gd3+, Er3+, Pb 6p² en acoplamiento jj con datos NIST, huecos, desdoblamiento
en Oh, magnetismo); efecto Zeeman; orbitales hidrogenoides y operadores; entrenador
«¿qué orbital es esta función?»; ensayo a la llama (práctica); Beer–Lambert
(práctica 4); color, fluorescencia y energía de fotones; radiactividad; RMN básica.
Diagrama de Grotrian según la etapa 9.
```

### Etapa 4 — Moléculas y enlace

```text
Vista «Moléculas y enlace» con un laboratorio unificado de ~60 moléculas: vistas
«Lewis + 3D» (por defecto), «3D grande», «Lewis plano», «Plegar» (animación de la
estructura plana a la 3D) y «Comparar» (dos moléculas lado a lado y por qué difieren);
cargas formales, hibridación, VSEPR, polaridad y momento dipolar, pares libres y
radicales (un electrón = un solo punto), grupo puntual y elementos de simetría,
energías de enlace, frecuencias de vibración experimentales; botón «vista frontal».
Etiquetas de geometría: negrita = experimental (con fuente), cursiva «≈» = modelo;
ángulos fijados por simetría (90, 109,47, 120, 180) se dan como exactos.
Galería 3D de solapamiento de orbitales (σ, π, δ, enlaces de hidrógeno, dímero de agua
O···O ≈ 298 pm, D0 = 13,2 kJ/mol) con nube electrónica, δ±, cargas formales, dipolo,
distancias y una ficha por escena. Diagramas de OM de diatómicas; Hückel para
aromáticos; unión de octaedros y tetraedros; boranos (reglas de Wade) y silanos.
Comprueba cada geometría con un análisis de grupo puntual.
```

### Etapa 5 — Cristalografía y estado sólido

```text
Vista «Cristalografía»: 7 sistemas y 14 redes de Bravais en 3D; índices de Miller
(planos y direcciones, entrenador); celda unidad (átomos por celda, fórmula, NC,
huecos, fracción de empaquetamiento, densidad); indexar un difractograma; deducir la
fórmula del cristal; conteo en redes 2D; calculadoras de los problemas de examen;
DRX y FRX. Vista «Estado sólido y materiales»: de OM a bandas, semiconductores y
dopaje, Ellingham, metalurgia y minerales de Bolivia, lantánidos y actínidos,
perovskitas (factor de tolerancia), baterías de ion litio, zeolitas/MOF/BET,
técnicas de caracterización. Visor 3D de estructuras (NaCl, CsCl, ZnS blenda y
wurtzita, CaF2, TiO2 rutilo, perovskita, espinela, corindón, diamante, grafito…),
boranos, silanos y moléculas. Poliedros. Ruta atrano y materiales porosos (MCM-41,
TEM sintética con transformada de Fourier, Tauc).
```

### Etapa 6 — Inorgánica I y química general

```text
Vistas: «Ácido–base, redox y grupos» (Arrhenius → Lewis, HSAB de Pearson, fuerza de
oxoácidos/hidruros/acuoiones, Born–Haber, Born–Landé y Kapustinskii, Latimer, Frost,
Pourbaix, química descriptiva de los grupos principales, defectos);
«Redox, termo y periódica» (igualador de ecuaciones por ion-electrón con análisis
redox, termodinámica, potenciales, propiedades periódicas con gráficos).
```

### Etapa 7 — Coordinación, organometálica y bioinorgánica (Inorgánica III)

```text
Vista «Coordinación y organometálica»: nomenclatura IUPAC 2005 de complejos (generador
con explicación); isomería con enumeración automática y modelos 3D; constantes de
formación, efecto quelato y diagramas de especiación; efecto trans y síntesis dirigida
de isómeros de Pt(II); OM de complejos octaédricos; regla de los 18 electrones;
carbonilos (retrodonación π e IR de ν(CO)); ciclos catalíticos (Wilkinson,
Monsanto/Cativa, Heck, hidroformilación); bioinorgánica (hemoglobina, citocromos,
clorofila, cobalamina, nitrogenasa, cisplatino).
```

### Etapa 8 — Técnicas instrumentales

```text
Una vista por técnica, cada una con: qué ocurre en el átomo/molécula, cómo es el
equipo (esquema), simulador del espectro, interferencias, y «Desde el curso:
problemas de examen»: Rayos X (DRX y FRX: espectro del tubo, Moseley, Auger, Compton,
Bragg, anillos de polvo y red recíproca, términos de las líneas K/L); Espectroscopía
atómica (AA, emisión, ICP, fluorescencia; ensanchamiento Doppler y Lorentz;
población de Boltzmann; elección de la línea de trabajo con datos de sensibilidad);
UV-Vis y fluorescencia (Jablonski, espectros de complejos); IR y Raman (modos normales
animados, oscilador armónico, efecto isotópico); RMN y EPR; Masas y GC-MS (sector
magnético, patrones isotópicos, fragmentación); Nuclear, gamma y Mössbauer del ⁵⁷Fe;
XPS y UPS; Microscopía SEM y TEM (longitud de onda del electrón, volumen de
interacción, SAED); Procesos industriales (Faraday, electrólisis).
```

### Etapa 9 — Ion libre frente a complejo: transiciones de s a g y color de d⁰ a d¹⁰

```text
(a) En «Átomo y espectros», diagrama de Grotrian para H, Li, Na, K, Cs, Cu, Ag, He, Mg,
Ca, Zn, Cd y Hg con energías del NIST Handbook. Columnas s, p, d, f, g (singletes y
tripletes para los de dos electrones; columna aparte para los niveles con el núcleo
excitado de Cu 3d⁹4s² y Ag 4d⁹5s²). Niveles no medidos se estiman con el defecto
cuántico (Rydberg–Ritz) y se dibujan discontinuos con «≈». Mostrar TODAS las
transiciones con colores distintos: permitidas (color de su λ), Laporte/misma
paridad (rosa), Δl o ΔL > 1 (naranja), intercombinación ΔS ≠ 0 (verde, fuerte en
Hg y Cd), J 0↔0 o ΔJ > 1 (gris); casillas para mostrar cada tipo; clic en un nivel →
todas sus transiciones; tabla con λ en aire (corrección de Edlén), energía, región;
modos absorción/emisión; población de Boltzmann.
(b) En «Campo cristalino», sección «¿Por qué hay color? De d⁰ a d¹⁰»: selector d⁰…d¹⁰
(espín alto/bajo en d⁴–d⁷), 5 orbitales degenerados del ion libre frente a t2g/eg con
los electrones, flecha del salto (permitido / solo invirtiendo espín / imposible),
espectro con franja visible, ion de ejemplo, color, bandas y ε (Lever 1984):
d⁰ Sc³⁺ incoloro (contraejemplo MnO4− violeta por TC); d¹ [Ti(H2O)6]3+ 20 300 cm⁻¹;
d² [V(H2O)6]3+ 17 800/25 700; d³ [Cr(H2O)6]3+ 17 400/24 600; d⁴ [Cr(H2O)6]2+ ~14 000;
d⁵ [Mn(H2O)6]2+ rosa pálido, ε ≈ 0,01–0,05; d⁶ [Fe(H2O)6]2+ 10 400 (bajo espín
[Co(NH3)6]3+ 21 000/29 500); d⁷ [Co(H2O)6]2+ 8 100/19 400/21 600; d⁸ [Ni(H2O)6]2+
8 500/13 800/25 300; d⁹ [Cu(H2O)6]2+ ~12 500; d¹⁰ Zn²⁺, Cd²⁺, Cu⁺, Ag⁺ incoloros
(contraejemplos AgI, Cu2O: banda prohibida del sólido o TC). Enlazar ambas secciones.
```

### Etapa 10 — Exámenes, entrenamiento y herramientas

```text
Vista «Exámenes y prácticas» con los datos del Anexo A: mapa de temas con frecuencia,
problemas resueltos, simulacro de examen con tiempo, avance por tema (guardado local),
hoja de fórmulas imprimible, preguntas por examen enlazadas a la sección de la app,
prácticas QCI 224 (Anexo B) con su laboratorio. Vista «Entrenamiento por niveles»
(preguntas con retroalimentación, «desmonta las respuestas rutinarias», estimaciones
de Fermi). Vistas «Guía y catálogo», «Símbolos y notación» y «Glosario».
Vista «Temas emergentes»: MOF (Nobel 2025), litio de salmueras (Uyuni), arsénico en
aguas, hidrógeno y electrocatálisis, puntos cuánticos, química verde, tierras raras,
baterías más allá del litio. Para temas posteriores a tu fecha de corte, busca en la
web y cita.
```

---

## 4. Acotaciones y lecciones aprendidas (sirven para futuras aplicaciones)

### Sobre el proceso de trabajo
1. **Construir por etapas y verificar cada una.** Cada versión se publicó, se adjuntó y se probó. Pedir «todo de una vez» produce páginas bonitas pero con contenido delgado.
2. **Siempre adjuntar el archivo además de publicar.** En la versión 1 solo se publicó y no se adjuntó, y el usuario no encontró el archivo.
3. **Distinguir «analiza» de «haz».** Cuando se pide analizar primero, responder con qué se perdió, por qué y qué se propone, sin tocar el código. Esto evitó cambios que el usuario no quería.
4. **No perder contenido al reorganizar.** Al rehacer vistas, algunas secciones (el diagrama de Grotrian, la galería) desaparecieron de su lugar original. La solución fue inventariar lo que se movió, rescatarlo con criterio químico y dejar enlaces desde el lugar antiguo.
5. **Rescatar de versiones antiguas con criterio.** De versiones anteriores se recuperó lo que era químicamente correcto y se descartó lo que no lo era (geometrías inventadas, nombres no IUPAC, animaciones sin base física), explicando cada decisión.
6. **Ser honesto con lo no verificado.** Cada entrega termina diciendo qué datos vienen de memoria y deben verificarse. Ejemplos: la interpretación de [Cu(CF₃)₄]⁻ sigue en debate; d⁷ de espín bajo no tiene un ejemplo octaédrico típico.

### Sobre «sin conexión»
7. **Three.js por CDN falla sin internet.** La primera versión del laboratorio 3D dependía de un CDN. La solución definitiva fue un **motor 3D propio en Canvas 2D** (rotación, zoom, esferas sombreadas, enlaces, profundidad), que funciona en cualquier navegador sin WebGL. WebGL queda como mejora opcional.
8. **Las bibliotecas útiles se guardaron en una «Biblioteca web de química»** (artefacto almacén) para incrustarlas en futuros proyectos: 3Dmol.js, OpenChemLib, RDKit (JS + WASM), KaTeX, SmilesDrawer y Three.js con OrbitControls. Las licencias son permisivas (BSD/MIT). Incrustarlas aumenta mucho el tamaño, así que solo se hace cuando aportan algo que el código propio no da.
9. **GPL/AGPL: no copiar.** Por ejemplo, Jmol/JSmol (LGPL) y otros visores: se estudió la idea y se reimplementó.

### Sobre la química (errores reales que se corrigieron)
10. **H₂O₂:** el ángulo diedro estaba en −10°. Al comprobarlo con pymatgen se vio el error, y se corrigió a 111,5° (fase gas).
11. **Grupos puntuales:** se comprobaron las 60 moléculas con PointGroupAnalyzer. N₂F₄, N₂O₅ y H₃PO₄ tienen geometría idealizada y se marcan así en la app.
12. **Radicales (NO₂):** el electrón desapareado se dibujaba como un par. Ahora es un único punto.
13. **Retos de geometría:** el reto mostraba la geometría real antes de responder. Ahora parte de la estructura plana canónica.
14. **Valores experimentales frente a modelo:** las etiquetas no los distinguían. Ahora la negrita es experimental, con fuente, y la cursiva «≈» es modelo.
15. **Nomenclatura:** «eteno» es el nombre preferido IUPAC 2013; «etileno» va entre paréntesis. «Acetileno» es un nombre retenido.
16. **Ion libre frente a complejo:** el diagrama de Grotrian es del átomo libre, donde los d están degenerados y d → d está prohibido por Laporte. Por eso no muestra desdoblamiento de campo cristalino; eso se ve en la sección d⁰–d¹⁰. Se aclararon cuatro confusiones frecuentes:
    - d¹ sí tiene color.
    - d⁰ y d¹⁰ no tienen bandas d–d.
    - d⁵ de espín alto es casi incoloro.
    - Cu²⁺ es d⁹, mientras que Cu⁺ y Ag⁺ son d¹⁰.
17. **Términos:** la expresión regular que leía los términos fallaba con letras posteriores a K (términos f). Se corrigió.
18. **Intercombinación:** en átomos pesados (Hg 253,7 nm, Cd 326,1 nm) es intensa. No se debe dibujar como «prohibida» igual que en He.

### Sobre el código
19. **Un solo IIFE grande:** conviene vigilar los nombres duplicados de nivel superior (hubo un choque `symElems` entre dos laboratorios). Antes de fusionar, buscar duplicados.
20. **Editar el bloque completo entre marcadores,** no fragmentos sueltos. Los reemplazos parciales dejaban versiones desincronizadas.
21. **`[hidden]{display:none!important}`:** sin esta regla, los elementos con `display:flex` ignoran `hidden`.
22. **Móvil:** las fichas y tablas anchas desbordaban (405 px en una pantalla de 390 px). Se resolvió con `table-layout:fixed`, `overflow-wrap:anywhere` y `.tablewrap{overflow-x:auto}`, y se añadió una prueba automática de desborde por vista.
23. **Leyendas SVG:** calcular la posición según el largo del texto o usar dos filas. Con un espaciado fijo los rótulos se pisaban.
24. **JavaScript:** `-x**2` es un error de sintaxis; hay que escribir `-(x**2)`.
25. **Rótulos SVG largos** se cortan en el borde del `viewBox`: dejar margen o partirlos en dos líneas.
26. **Enlaces entre vistas:** el router debe aceptar ids de secciones de otras vistas, cambiar de vista y luego desplazarse, con un pequeño retraso porque la vista se renderiza de forma perezosa.
27. **Publicación con conflicto de versión:** leer la versión en vivo, compararla y solo entonces publicar sobre ella.

### Sobre la enseñanza
28. **Cada laboratorio responde tres preguntas:** qué pasa, por qué pasa y qué pasaría si cambio X.
29. **Mostrar también lo prohibido y lo que no ocurre.** Las transiciones prohibidas en otro color y los casos sin color (d⁰, d¹⁰) enseñan tanto como los permitidos.
30. **Conectar teoría, examen y práctica:** cada tema enlaza a sus problemas de examen y a su práctica de laboratorio, y cada pregunta de examen enlaza a la sección donde se estudia.
31. **Contexto local:** minerales y metalurgia de Bolivia, litio del salar de Uyuni y arsénico en aguas del altiplano.

---

### Sobre la versión «Laboratorio Interactivo de Química Inorgánica 3D» (análisis de un prompt externo)
32. **De un prompt externo se rescata el método, no el contenido.** El prompt de mejora que se analizó era de otra aplicación (ShEPhERD-score, curvas ROC); se tomaron sus principios de diseño y pedagogía y se descartó el resto.
33. **Dos tablas de colores de elementos** convivían con tonos distintos para C y N. Se unificó en una sola tabla Jmol/CPK.
34. **Una fuente nueva no puede romper el modo sin conexión:** se incrustan en base64 (unos 150 KB en total). Lexend no tiene griego ni subíndices Unicode; el navegador los toma de la fuente del sistema.
35. **Quitar las barras fijas empeora la navegación** en una aplicación de más de 30 vistas: la solución intermedia es ocultarlas al bajar.
36. **Las afirmaciones de neurociencia se formulan con prudencia:** que predecir antes de observar ayuda está bien documentado; que una tipografía concreta mejore la lectura es más discutido.
37. **No obligar a todos a predecir:** bloquear cada laboratorio estorba a docentes y evaluadores. Lo mejor es que cada usuario elija el modo, con «Explorar» por defecto.
38. **Antes de añadir un desplegable para principiantes**, comprueba si la sección ya tiene uno, para no duplicarlo.
39. **Al citar caracteres de una tabla de grupos**, comprobar el valor: el de los orbitales d bajo C₄ en O_h es −1, no +1.
40. **Comprueba las constantes al cambiar de unidades:** la constante de Kapustinskii es 120,25 kJ nm mol⁻¹ con distancias en nm, o 1213,8 con distancias en Å.
41. **Datos sin servidor:** el reporte lo genera y envía el propio estudiante. Así se mantiene el funcionamiento sin internet y no se guardan datos personales en ningún sitio. La firma solo detecta ediciones a mano: no es seguridad fuerte.
42. **Grado de seguridad junto a cada respuesta:** un error con «seguro» revela una idea previa equivocada, y es lo que más le sirve al docente.
43. **En la página publicada en claude.ai, el visor bloquea las descargas:** por eso conviene ofrecer también «copiar el texto» e imprimir.

44. **Créditos neutrales:** no se debe atribuir la obra a una cátedra o a un docente que no la avala. Basta con el autor, las fuentes (el plan de estudios y la bibliografía) y un agradecimiento sobrio, con la aclaración de que no es una publicación oficial. Así nadie se siente comprometido a participar o a no hacerlo.

## Anexo A — Exámenes pasados de Química Inorgánica I (sustituye al PDF)

**28 exámenes:** 2015-11-30, 2019-09-18, 2019-10-28, 2019-12-16, 2020-06-05, 2020-07-17, 2020-07-24, 2020-07-27 (×3), 2020-08-14, 2020-10-07, 2020-11-09, 2020-12-09, 2020-12-16, 2021-04-01, 2021-05-13, 2021-06-10, 2021-06-15, 2021-09-27, 2021-11-01, 2021-12-06, 2022-04-13, 2022-05-18, 2022-06-15, 2022-06-20, 2022-09-26 y 2022-11-07.

**Frecuencia de los temas** (número de preguntas):

| Tema | Veces | Sección de la app |
|---|---|---|
| Celda unidad: densidad, Z, radios | 35 | Cristalografía → calculadoras |
| Lewis, cargas formales, hibridación, geometría, polaridad | 27 | Moléculas |
| Difracción de rayos X | 18 | Cristalografía / Rayos X |
| RMN | 18 | Átomo / RMN y EPR |
| Términos y microestados | 16 | Átomo |
| Fórmula a partir de NC y conteo de átomos | 16 | Cristalografía 2D |
| Fluorescencia de rayos X | 14 | Rayos X |
| Índices de Miller | 14 | Cristalografía |
| Recuento de OM σ y π | 12 | Moléculas |
| Espectrometría de masas | 11 | Masas |
| Radiactividad y núcleo | 11 | Nuclear |
| Redox y procesos industriales | 10 | Redox / Industriales |
| Termodinámica | 10 | Redox y termo |
| Boranos | 10 | Moléculas |
| Huecos octaédricos y tetraédricos | 10 | Cristalografía |
| Espectro del hidrógeno | 9 | Átomo |
| Campo cristalino | 9 | Campo cristalino |
| Orbitales y operadores | 9 | Átomo |
| Sistemas cristalinos | 6 | Cristalografía |
| Color, UV-Vis y fotones | 5 | Átomo / UV-Vis |
| Diagramas de OM de diatómicas | 5 | Moléculas |
| Silicatos y polímeros inorgánicos | 5 | Cristalografía |
| Sodio y efecto Zeeman | 4 | Átomo |
| Propiedades periódicas | 3 | Periódica |

**Ejemplo de enunciado típico** (tercer parcial, 2015):
- Electrólisis de salmuera y de alúmina: ecuaciones, oxidante y reductor.
- FRX del cobre: asignar 8,028, 8,048, 8,905 y 0,930 keV a Kα1, Kα2, Kβ y Lα; dar el estado inicial y final, y la λ de Kα.

## Anexo B — Prácticas de laboratorio QCI 224 (sustituye al .docx)

1. **Osciloscopio:** VOLTS/DIV, TIME/DIV, disparo, acoplamiento AC/GND/DC, modo X-Y, calibración, voltaje de pilas, ondas sinusoidales y cuadradas (25–200 Hz), constante de tiempo RC (R = 33 000 Ω, C = 0,049 μF).
2. **Figuras de Lissajous:** relaciones 1:1 y 1:2; x = A sen(ωt), y = B sen(nωt + φ).
3. **Absorción y emisión atómicas:** líneas del hidrógeno (Rydberg), doblete del sodio y efecto Zeeman.
4. **Absorción molecular:** ley de Beer, color complementario, transiciones d–d (Tanabe–Sugano) y π→π*.
5. **Vida media del tritio y del ¹⁴C:** A = λN, centelleo líquido (CPM y eficiencia), ln(A₀/A) = λt.
6. **RMN:** espín nuclear, frecuencia de Larmor ν = (γ/2π)B, desplazamiento químico, multiplicidad n + 1.
7. **Estructuras compactas y no compactas:** sc, bcc, fcc (ccp) y hcp; NC; huecos.
8. **Volúmenes de ocupación:** sc 52,4 %, bcc 68,0 %, fcc y hcp 74,0 %, diamante 34,0 %.
9. **Relación de radios:** límites 0,155 (NC 3), 0,225 (4), 0,414 (6), 0,732 (8) y 1 (12).
10. **Sistemas cristalinos 2D y 3D:** conteo en redes 2D (1/4, 1/2, 1); 7 sistemas y 14 redes de Bravais.
11. **Fluorescencia de rayos X:** líneas K y L, ley de Moseley, cristal analizador y 2θ.
12. **Difracción de rayos X:** ley de Bragg, índices de Miller, extinciones sistemáticas, densidad a partir de la celda.

---

*Página publicada de referencia:* https://claude.ai/artifact/TwbpHav4RSczpFPYX7Ni7K
