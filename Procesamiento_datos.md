# Procesado de datos metagenómicos
El análisis bioinformático de datos metagenómicos tiene como objetivo estudiar la diversidad microbiana y las funciones potenciales presentes en una muestra ambiental (como suelo, agua, microbioma humano, etc.) sin necesidad de cultivar los microorganismos. Este tipo de análisis se basa en secuenciación masiva (NGS) y sigue una serie de pasos bien definidos que garantizan la calidad, precisión y utilidad biológica de los resultados obtenidos.
El flujo de trabajo modular (calidad → filtrado de hospedero → deduplicación → QC → ensamblado/clasificación → perfil funcional → análisis de virulencia) es la espina dorsal de la mayoría de pipelines metagenómicos shotgun y puede adaptarse según volumen de datos, objetivos del proyecto y recursos computacionales disponibles.

## 1. Control de calidad y recorte de adaptadores
El análisis comienza con la eliminación de secuencias de baja calidad y de los adaptadores o primers específicos que se hayan añadido durante la construcción de la librería. Herramientas como Trimmomatic o Trim Galore se emplean con parámetros que definen el umbral de calidad de los nucleótidos (usualmente Phred ≥ 33), una longitud mínima de lectura (por ejemplo, 100 bp) y una “stringencia” en la coincidencia del adaptador. Este paso asegura que sólo las lecturas que superan la calidad mínima entran en las etapas posteriores, reduciendo el ruido y mejorando la fiabilidad del ensamblado o de los perfiles taxonómicos y funcionales.
```Bash
# FastQC (multithreaded)
fastqc -t 8 sample_R1.fastq sample_R2.fastq -o qc_reports/

# Trimmomatic (paired-end, multithreaded)
trimmomatic PE -threads 8 \
  sample_R1.fastq sample_R2.fastq \
  sample_R1.trimmed.fastq sample_R1.unpaired.fastq \
  sample_R2.trimmed.fastq sample_R2.unpaired.fastq \
  SLIDINGWINDOW:4:20 MINLEN:100

# Trim Galore (paired‐end, multithreaded)
trim_galore --paired --cores 8 \
  --length 100 --stringency 3 \
  sample_R1.fastq sample_R2.fastq \
  -o trimmed/
```
## 2. Eliminación de secuencias del hospedero
Cuando la muestra proviene de un organismo anfitrión (por ejemplo muestras de murciélagos o humanas), es fundamental filtrar esas lecturas antes del análisis microbiano. Para ello, se descargan las referencias genómicas del hospedero (p. ej. Myotis myotis, R. ferrumequinum) y se alinean las lecturas contra ellas usando Bowtie2. Las lecturas que mapean al genoma hospedador se descartan, evitando que el ADN “host” contamine los resultados finales y sesgue la composición microbiana y funcional.
También se pueden eliminar contaminación ambiental o artefactos utilizando bases de datos personalizadas o filtros manuales.

El sitio web https://www.ncbi.nlm.nih.gov/datasets/genome/ forma parte de los servicios ofrecidos por el NCBI (National Center for Biotechnology Information) y está diseñado para facilitar el acceso a datos genómicos de organismos de interés en investigación biomédica, agrícola o ecológica. 
Esta plataforma permite a los usuarios buscar, descargar y analizar secuencias completas de genomas, así como metadatos asociados, desde una interfaz intuitiva y estandarizada. Además, ofrece herramientas como el Datasets API y el CLI (Command Line Interface) , que permiten la descarga automatizada de archivos en formato FASTA, GenBank, GFF y otros, lo cual es especialmente útil en análisis comparativos, reconstrucción filogenética o estudios funcionales. 
El servicio también se integra con otras bases del NCBI, como BioProject, SRA y Protein, convirtiéndose en un recurso centralizado para el estudio de la genómica funcional y evolutiva.
Desde esta plataforma se pueden extraer los genomas de interés para construir los índices necesarios durante el procesamiento de datos.

```Bash
# Index build
bowtie2-build --large-index Mixed.fasta host_genome_index

# Bowtie2 (host‐read removal, multithreaded)
bowtie2 -x host_genome_index -1 sample_R1.trimmed.fastq -2 sample_R2.trimmed.fastq \
  -p 8 --very-sensitive -S host_mapped.sam --un-conc-gz host_unmapped \
  && samtools view -b -f12 host_mapped.sam > host_removed.bam
```
## 3. Eliminación de duplicados y artefactos
Muchas veces, por la PCR o por errores de flujo de trabajo, aparecen lecturas duplicadas o artefactos como homopolímeros extremos (por ejemplo cadenas largas de “G”). Con utilidades como PRINSEQ-lite se detectan y eliminan duplicados (--derep) y lecturas que no cumplen un mínimo de calidad media. Además, Trim Galore u otras herramientas pueden eliminar regiones de homopolímeros excesivos mediante un adaptador artificial, evitando que estos fragmentos afecten los posteriores ensamblados o cuantificaciones.
```Bash
# PRINSEQ‑lite (deduplicación; no multithreading)
perl prinseq-lite.pl -fastq host_removed.fastq \
  -derep 1 -min_len 100 -min_qual_mean 20 -out_format 3 \
  -out_good deduped -out_bad discarded
```
## 4. Control de calidad previo y posterior
Tanto al inicio como después de cada gran filtrado, se corre un informe de calidad con FastQC. Este paso sirve para monitorizar métricas clave (distribución de calidad por posición, contenido de GC, presencia de adaptadores remanentes) y verificar que las operaciones de limpieza están funcionando correctamente. Un aumento notable de lecturas cortas, un repunte de homopolímeros o un desequilibrio en la distribución de bases pueden indicar problemas que requieren ajuste de parámetros.
```Bash
# Data base download
wget https://genome-idx.s3.amazonaws.com/centrifuge/p%2Bh%2Bv.tar.gz
tar -xvzf p+h+v.tar.gz

# Centrifuge (classification, multithreaded)
centrifuge -x p+h+v -1 deduped_1.fastq -2 deduped_2.fastq \
  -p 8 --report-file centrifuge.report \
  -S centrifuge.out
```
## 5. Ensamblado o clasificación taxonómica
Con las lecturas filtradas, existen dos rutas principales:
Ensamblado (de novo) con programas como MEGAHIT o MetaSPAdes, para reconstruir genomas completos o contigs representativos.
Clasificación directa con herramientas basadas en k‑mers (p. ej. Centrifuge) que asignan cada lectura a un taxón usando bases de datos preconstruidas (bacterias, virus, hongos). Esta segunda opción es más rápida y permite cuantificar la abundancia relativa de cada organismo.

```Bash


```
## 6. Perfil funcional y análisis de rutas metabólicas
Más allá de “quién está” en la microbiota, interesa saber “qué hacen”. Para ello, se usan suites como HUMAnN, que mapean las lecturas contra colecciones de genes de referencia y anotan rutas metabólicas (p. ej. vía UniRef y MetaCyc). El resultado es un perfil de abundancia de funciones y rutas, muy útil para comparar ecologías microbianas entre condiciones experimentales.
```Bash
# HUMAnN 3.0 (functional profiling, multithreaded)
humann --input deduped.fastq --output humann_out/ --threads 8

```

## 7. Detección de factores de virulencia
Para estudios de patogenicidad o vigilancia epidemiológica, se buscan genes relacionados con virulencia usando BLAST contra bases de datos especializadas (por ejemplo VFDB). Con umbrales estrictos de identidad (≥ 95 %) y valores de e‑value (≤ 1e‑10), se identifican potenciales factores de virulencia presentes en la muestra.
```Bash
# Download data base
wget https://ftp.ncbi.nlm.nih.gov/blast/db/VFDB.tar.gz
tar -xvzf VFDB.tar.gz
 
# Index build
makeblastdb -in VFDB.fasta -dbtype nucl -out VFDB

# BLASTn (virulence search, multithreaded)
blastn -query deduped.fastq -db VFDB \
  -perc_identity 95 -evalue 1e-10 \
  -num_threads 8 -out vf_hits.tsv
```
## Datos que debemos registrar de acuerdo a formulario MIxS
Las secuencias genómicas deben formatearse de acuerdo con el Consorcio de Estándares Genómicos (GSC) y seguir la especificación de información mínima sobre una secuencia x (MIxS). 

MIxS "Minimal Information about (X) Sequence": Es un conjunto de estándares desarrollados por el Genomic Standards Consortium (GSC) para describir datos de secuenciación con metadatos mínimos obligatorios y opcionales. La X puede ser:

* Genoma (MIGS),
* Metagenoma (MIMS),
* Genoma ambiental (MIENS),
* u otras variantes.

Garantiza que los datos de secuenciación vengan acompañados de metadatos relevantes, como lugar y fecha de muestreo, condiciones ambientales, plataforma de secuenciación, tipo de muestra (agua, suelo, humano, etc.). Facilita la reproducibilidad, comparación entre estudios y el almacenamiento estandarizado en bases de datos públicas como ENA, NCBI o MGnify.

Para más información ver:
https://www.gensc.org
