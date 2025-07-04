# Obtención de datos y gestión de dependencias

## Conda para el Análisis Bioinformático
En el ámbito del análisis bioinformático, la ejecución de herramientas especializadas suele requerir múltiples dependencias: bibliotecas específicas, versiones controladas de programas, entornos compatibles y configuraciones precisas. Este escenario puede volverse complejo rápidamente si se trabaja con múltiples proyectos que requieren distintas versiones de las mismas herramientas.

Conda surge como una solución eficiente y robusta para gestionar este tipo de desafíos. Se trata de un gestor de paquetes y entornos multiplataforma (compatible con Windows, Linux y macOS), capaz de instalar, administrar y aislar software en entornos independientes, garantizando así la reproducibilidad y estabilidad de los análisis científicos.

Una de las principales ventajas de Conda es su capacidad para manejar no solo paquetes de Python, sino también herramientas escritas en otros lenguajes como R, Perl, C++, entre otros. Esto lo convierte en una herramienta fundamental en bioinformática, donde muchas aplicaciones críticas (como samtools, bwa, hisat2, trimmomatic, prokka, etc.) están basadas en binarios compilados o scripts de sistemas Unix.

Además, Conda permite crear entornos virtuales personalizados, donde cada proyecto puede tener sus propias versiones de software sin interferir con otros. Esto es clave en contextos científicos, donde el uso de versiones específicas de herramientas puede ser crucial para asegurar la validez de los resultados.

Gracias a su integración con repositorios como Bioconda, un canal especializado en software bioinformático, Conda ofrece acceso directo a miles de herramientas ya empaquetadas y listas para usar, facilitando enormemente la implementación de pipelines de análisis genómico, transcriptómico, metagenómico y más.

Antes de continuar, nos aseguraremos de tener herramientas básicas que nos servirán en la monitorización de nuestros recursos:

```Bash
# Incluye compiladores y herramientas necesarias para compilar codigo desde fuentes
sudo apt install build-essential

# Herramientas para descargar archivos desde la web.
sudo apt install curl wget

# Editores de texto para la linea de comandos.
sudo apt install vim nano

# Herramienta para monitorizar el uso de CPU, memoria y procesos en tiempo real.
sudo apt install htop
sudo apt install glances

# Utilidades para comprimir y descomprimir archivos.
sudo apt install zip unzip
sudo apt install pigz

# Algunos programas bioinformaticos requieren Java.
sudo apt install default-jdk

# Algunos scripts bioinformaticos pueden estar escritos en Perl.
sudo apt install perl
```

Anaconda, Miniconda, y Bioconda son herramientas relacionadas con el ecosistema de Conda, pero tienen funciones y alcances distintos. 
Anaconda es una distribución completa que incluye Conda, Python, y cientos de paquetes científicos preinstalados, orientada a facilitar el trabajo en ciencia de datos, machine learning e investigación. 
Debido a su tamaño, suele ser más pesada de instalar. Miniconda, en cambio, es una versión ligera que solo incluye Conda y Python, permitiendo al usuario instalar únicamente los paquetes que necesite, lo cual lo hace ideal para entornos controlados como la bioinformática. 
Finalmente, existen diferentes versiones de Conda, algunas enfocadas en rendimiento (como Mamba), y otras adaptadas a necesidades específicas, lo que permite flexibilidad y escalabilidad en proyectos científicos.


```Bash
# Instalar Anaconda
cd ~
wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
bash Anaconda3-2024.06-1-Linux-x86_64.sh # ejecuta el instalador
source ~/.bashrc # reinicia la terminal

# verifica instalación
conda --version
python --version

# Instalar Miniconda
cd ~
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
bash Miniconda3-latest-Linux-x86_64.sh # ejecuta el instalador
source ~/.bashrc # reinicia la terminal

# verifica instalación
conda --version
python --version
```

## Canales, ambientes y compartimentación en Conda
### ¿Qué son y para qué sirven los canales en Conda?
Los canales (channels ) son fuentes oficiales o comunitarias donde se almacenan y distribuyen paquetes compatibles con Conda. 
Actúan como repositorios de software, permitiendo que los usuarios instalen herramientas desde diferentes fuentes, dependiendo del tipo de proyecto o necesidad.

Algunos canales populares:
* defaults: El canal oficial por defecto que viene con Anaconda o Miniconda. Incluye paquetes esenciales de Python, R y algunas herramientas científicas.
* conda-forge: Un canal comunitario altamente actualizado y amplio, mantenido por la comunidad de Conda. Ofrece versiones recientes de muchos paquetes y a menudo tiene mayor cobertura que defaults.
* bioconda: Un canal especializado en software bioinformático. Contiene miles de herramientas útiles para análisis genómico, transcriptómico, metagenómico, entre otros.
* pytorch, tensorflow, r, etc.: Canales dedicados a ecosistemas específicos como machine learning o estadística avanzada.

El orden en que se especifican los canales a la hora de configurar conda afecta la prioridad de descarga de los paquetes de interés. 

```Bash
# Configurar canales
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Buscará primero en bioconda y luego en conda-forge.

### ¿Qué son los ambientes?
Los ambientes (environments) son espacios aislados donde se instalan paquetes y sus dependencias sin interferir con otros proyectos. 
Esto permite manejar distintas versiones de un mismo programa o lenguaje según las necesidades de cada proyecto. Por ejemplo puedes tener un ambiente con R 4.1 y otro con R 4.3.

Los ambientes ofrecen una forma de aislar dependencias y versiones de software, lo cual es fundamental, donde diferentes proyectos pueden requerir versiones distintas de una misma herramienta.
También es idea cuando es importante garantizar la reproducibilidad de los análisis. La compartimentación evitan conflictos entre bibliotecas y dependencias cruzadas. 
Otra de las ventajas de la compartimentación es que se puede compartir un entorno exacto con otros investigadores mediante un archivo environment.yml.

```Bash
# Comando para crear ambiente conda llamado "metabiome-preprocessing":
conda create -n metabiome-preprocessing fastqc multiqc trimmomatic samtools bowtie2 trim-galore bamtools prinseq

# Crear el ambiente a partir de un archivo YAML:
conda env create -f metabiome-preprocessing.yml

# Actualizar todos los paquetes:
conda update --all

# Cuando termines de trabajar con el entorno, puedes desactivarlo:
conda deactivate

# Extraer el ambiente conda:
conda env export > metabiome-preprocessing.yml
```
---
## Secuenciación para obtener datos

### Preparación de librerías para Illumina (metagenómica shotgun)
#### Paso 1: Extracción y calidad del ADN :
Primero se extrae ADN total del entorno (agua, suelo, heces…), utilizando kits según el tipo de muestra (p.ej., PowerSoil Kit para suelos). 
Se verifica la cantidad (por técnicas fluorométricas) y la pureza (UV con espectrofotometría: ratios A260/280 y A260/230 para detectar impurezas).

#### Paso 2: Fragmentación del ADN
El ADN genómico se fragmenta en trozos de 300–600 bp mediante:
Shearing mecánico (sonicación) — mayor uniformidad, menos sesgo.
Fragmentación enzimática (tagmentación, trasposasas) — proceso más rápido pero puede introducir sesgo.

#### Paso 3: Reparación de extremos y A‑tailing
Tras fragmentar, los extremos quedan "irregulares":
Se reparan y se añade un nucleótido ‘A’ en el extremo 3′ para facilitar la ligación del adaptador.

#### Paso 4: Ligación de adaptadores
Se unen adaptadores tipo P5 y P7, que contienen:
Secuencias complementarias a los oligos del flow cell.
Índices (barcodes) para barcoding y multiplexación (p. ej. hasta 96 muestras).
Sitios de anillamiento para primers.

#### Paso 5: Limpieza y selección de tamaño
Se eliminan fragmentos no deseados (p. ej., dímeros de adaptador) usando perlas magnéticas (AMPure XP) o electroforesis. 
Esto mejora la calidad final de la librería.

#### Paso 6: (Opcional) Amplificación por PCR
Cuando la cantidad inicial de ADN es baja, se realiza una PCR para enriquecer las moléculas ligadas con adaptador. 
Reduce sesgo, pero exagerarla puede introducir duplicados; en kits modernos se tiende a minimizar esta etapa.

#### Paso 7: Cuantificación y control de calidad final
Se mide la concentración (p. ej. Qubit, qPCR) y calidad (Bioanalyzer o Tapestation) antes de cargar sobre el flow cell. 

---
### Secuenciación por síntesis (Sequencing by Synthesis – SBS)
Una vez lista la librería, se carga en el flow cell del equipo Illumina (MiSeq, NextSeq, NovaSeq). El flujo es:

#### a) Amplificación clonal (Bridge o ExAmp)
Los fragmentos se unen a oligos en la superficie y se amplifican en clusters todos idénticos:
Bridge amplification en flow cells tradicionales
Exclusion amplification (ExAmp) en flow cells con diseño ‘patterned’ 

#### b) Secuenciación cíclica
Cada ciclo de SBS añade una base marcada con fluorescencia y un terminador reversible:
Se incorpora una base (A, C, G o T).
Se toma una foto del patrón luminoso.
El bloque se quita y se repite el ciclo.
Se pueden generar lecturas de hasta 2×251 bp (lectura pareada) 

#### Aplicación a shotgun metagenómica
En metagenómica shotgun se secuencia todo el ADN en la muestra sin sesgos de amplificación de regiones específicas. Permite:
Catalogar diversidad microbiana (bacterias, virus, hongos, etc.)
Detectar genes funcionales, resistencia antimicrobiana, metabolitos 
Requiere mayor profundidad de secuenciación (p. ej. cientos de millones de lecturas por muestra) que mapea la estructura funcional y taxonómica de comunidades complejas.

## Extracción de datos del Sequence Read Archive del National Center for Biotechnology Information

El Sequence Read Archive (SRA) es una base de datos pública administrada por el NCBI que almacena secuencias de ADN y ARN de alto rendimiento generadas por tecnologías como Illumina, PacBio o Nanopore.
Este repositorio se ha convertido en una fuente fundamental para la investigación genómica, permitiendo a científicos de todo el mundo acceder a datos experimentales de secuenciación masiva. 
Sin embargo, los archivos almacenados en formato binario comprimido .sra no son directamente utilizables para análisis, por lo que se requiere de herramientas especializadas para su descarga, conversión y procesamiento. 
La correcta extracción de estos datos permite acceder a las lecturas originales (en formatos como FASTQ) para realizar análisis posteriores como ensamblaje, alineamiento, cuantificación génica o detección de variantes.

Para extraer los datos del SRA, se utiliza comúnmente el SRA Toolkit, un conjunto de herramientas proporcionado por el NCBI diseñado específicamente para interactuar con los archivos .sra. 
Una vez instalado, comandos como fastq-dump permiten convertir los archivos .sra a formatos estándar como FASTQ, ampliamente usados en pipelines bioinformáticos. 
Además, existen alternativas más rápidas como fasterq-dump, parte del toolkit, que mejora el tiempo de procesamiento gracias a una arquitectura optimizada. 
Estas herramientas permiten también filtrar lecturas técnicas, limitar la cantidad de datos a extraer o incluso trabajar directamente desde el acceso web al repositorio, sin necesidad de descargar previamente el archivo completo. 

```Bash
# Instala el SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/Downloads
# Lo instalaremos desde conda
conda install -c bioconda sra-tools

# Activa el ambiente
conda activate sra-tools

```
El SRA organiza la información genómica en una estructura en forma de árbol, con varios niveles de anidamiento. 
Esta organización permite almacenar, acceder y compartir datos desde el nivel más general hasta el más específico, facilitando la navegación y reutilización de datos.

```Bash
Study (Estudio)
└── Sample (Muestra biológica)
    └── Experiment (Experimento tecnológico)
        └── Run (Ejecución de secuenciación)
```
**Study (Estudio / Proyecto):** Es el nivel más alto de la estructura. Representa un proyecto o estudio científico completo, como un artículo publicado o un experimento amplio.
Contiene información general como objetivo del estudio, autores, instituciones, etc. Puede contener múltiples muestras biológicas y condiciones experimentales

**Sample (Muestra biológica):** Representa una muestra biológica específica utilizada en el estudio. Se identifica mediante un acceso único.
Incluye metadatos como especie, tejido, condición, sexo, edad, etc.

**Experiment (Experimento tecnológico):** Describe cómo se preparó y secuenció una muestra particular. Incluye detalles sobre la biblioteca de secuenciación (tipo: RNA-seq, ChIP-seq, etc.),
Tipo de plataforma usada (Illumina, PacBio, Nanopore…), Estrategia de secuenciación (paired-end, single-end).

**Run (Ejecución de secuenciación):** Es el nivel más bajo y contiene los datos brutos de secuenciación. Cada ejecución corresponde a un archivo .sra que contiene millones de lecturas (reads).
Puede haber varias ejecuciones para un mismo experimento si se secuencia más de una vez o en diferentes flujos celulares (lanes).

*Ejemplo práctico:*

* Study: PRJNA123456 – “Microbioma en pacientes con enfermedad inflamatoria intestinal”
* Sample: SRS456789 – “Muestra fecal de paciente con colitis ulcerosa”
* Experiment: SRX789012 – “Secuenciación 16S rRNA con primers V3-V4 en Illumina MiSeq”
* Run: SRR987654 – “Datos brutos obtenidos en corrida de secuenciación”

```Bash
# Obten las Run Accessions desde la pagina de resultados del SRA, selecciona las muestras de interes.
# https://www.ncbi.nlm.nih.gov/sra/


# Utiliza el comando prefetch para descargar los archivos .sra:
prefetch SRR12345678

# Convierte los archivos .sra a formato FASTQ utilizando fasterq-dump:
fasterq-dump --threads 4 --split-files SRR...
```

## Notas técnicas sobre el rendimiento en el analisis bioinformático

Para procesar grandes volúmenes de datos metagenómicos, el procesador (CPU) es fundamental. Lo ideal es contar con al menos 8 núcleos físicos; si la CPU soporta hyperthreading, cada núcleo físico puede atender dos hilos lógicos, doblando la capacidad de procesamiento (por ejemplo, 16 hilos lógicos en un CPU de 8 núcleos). 
El multithreading es la capacidad de los programas de dividir su trabajo en varios hilos que corren a la vez, aprovechando al máximo todos esos núcleos e hilos lógicos.

La memoria RAM determina cuánto de tu proyecto puede “vivir” en memoria sin recurrir al disco. Para análisis rutinarios recomendamos un mínimo de 32 GB; para ensamblajes complejos o conjuntos de muestra muy grandes, 64–128 GB darán margen extra. 
La velocidad de la RAM (por ejemplo, DDR4 o DDR5 a 2400 MT/s o más) y el número de canales (dual, triple o quad–channel) influyen directamente en el ancho de banda de memoria y, por ende, en la rapidez con la que los núcleos acceden a los datos.

Cuando la RAM se agota, el sistema usa swap, un área de disco que hace las veces de memoria virtual. Aunque evita que los procesos mueran por falta de memoria, el acceso a swap es mucho más lento que a la RAM: de ahí que un exceso de uso de swap provoque caídas drásticas en el rendimiento.

El almacenamiento juega un rol clave: un SSD NVMe (con velocidades de lectura/escritura superiores a 2 GB/s) mantiene bajas las esperas de E/S y maximiza la utilización de la CPU. En cambio, un HDD SATA sólo conviene para archivos de respaldo o datos fríos, pues su acceso aleatorio lento dispara el iowait, es decir, el porcentaje de tiempo que la CPU permanece inactiva esperando al disco.

El iowait es un indicador directo de cuellos de botella de disco: valores por encima del 20 % sugieren que la CPU pasa demasiado tiempo a la espera de lectura/escritura. Reducir el iowait implica mejorar el subsistema de almacenamiento (más SSDs, RAID, o NVMe) o ajustar la planificación de tareas para no saturar la E/S.

Para la transferencia de datos (descarga de SRA, carga de resultados al cloud, etc.) se recomienda una conexión de, al menos, 100 Mbps simétricos y latencia baja (<50 ms). Esto agiliza la transferencia de muchos archivos pequeños y evita interrupciones en flujos de trabajo que dependen de datos remotos.

---

> **Requerimientos técnicos para análisis metagenómico**  
>  
> **CPU (procesador)**  
> - Cores (núcleos físicos): ≥ 8 (recomendado 16+)  
> - Hyperthreading (hilos lógicos): duplica los threads disponibles, útil para multitarea.  
> - Multithreading: capacidad de un programa para dividir tareas en múltiples hilos y aprovechar varios cores.  
>  
> **Memoria (RAM)**  
> - Mínimo 32 GB, preferible 64–128 GB para ensamblajes grandes.  
> - Tipo DDR4/DDR5: velocidad ≥ 2400 MT/s; canales dual/triple/quad-channel mejoran ancho de banda.  
> - Swap: espacio en disco que actúa como “memoria virtual” cuando la RAM se agota, pero **muy** más lento.  
>  
> **Almacenamiento (disco)**  
> - SSD NVMe: velocidades de lectura/escritura ≥ 2 GB/s, reduce I/O wait.  
> - HDD SATA: sólo para datos fríos o backup; lentos para acceso aleatorio.  
> - I/O wait: porcentaje de tiempo que la CPU espera a que el disco responda; valores altos (> 20 %) indican cuello de botella de disco.  
>  
> **Red (internet)**  
> - ≥ 100 Mbps simétricos si descargas datos (p. ej. SRA) / subes resultados al cloud.  
> - Latencia baja (< 50 ms) mejora eficiencia en transferencias de múltiples archivos pequeños.  
>  
> **Cuellos de botella comunes**  
> 1. **I/O**: disco lento → alta I/O wait → CPU ociosa.  
> 2. **RAM insuficiente**: exceso de swap → gran lentitud.  
> 3. **Threads mal dimensionados**: usar demasiados hilos en SSD o poca RAM degrada rendimiento.  
>  
> **Conceptos clave**  
> - **iowait**: tiempo que la CPU está inactiva esperando operaciones de E/S (I/O) de disco.  
> - **swap**: espacio en disco que extiende la RAM; evita errores por falta de memoria pero es muy lento.  
> - **multithreading**: ejecución concurrente de múltiples hilos en un programa para paralelizar tareas.  

Algunas herramientas que pueden ayudar a vigilar el rendimiento de nuestros procesos nos ayudan a mejorar la eficiencia en nuestros análisis. En nuestro caso utilizaremos dos:

```Bash
glances
htop
```

# **Ejercicio:** descargar y convertir a formato fastq 8 runs (aleatorios) del proyecto PRJNA1269778.
```Bash
prefetch SRR(de acuerdo a la elección)
fastq-dump --split-files (ajustar al número de nucleos del procesador) -X 10000 -O ./data/ SRR(de acuerdo a la elección)
```
