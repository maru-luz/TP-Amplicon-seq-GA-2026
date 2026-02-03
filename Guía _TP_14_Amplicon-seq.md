# Trabajo Práctico: Amplicon-seq para detección de isoformas con lecturas largas (ONT)

## 1. Introducción

El objetivo de este trabajo práctico es introducir el uso de **secuenciación de lecturas largas (Oxford Nanopore Technologies, ONT)** aplicada a **cDNA-amplicon-seq** para la detección y caracterización de **isoformas transcriptómicas**.

A diferencia de tecnologías de lecturas cortas (Illumina), las lecturas largas permiten observar directamente combinaciones completas de exones dentro de una misma molécula, lo que resulta clave para:

* Identificar isoformas completas
* Detectar eventos de splicing alternativo
* Distinguir isoformas canónicas de isoformas novedosas

Durante el TP se trabajará con datos reales de ONT correspondientes a amplicones del gen ***FMR1***, a partir de dos muestras biológicas (Sample A y Sample B) con un conjunto parcialmente compartido de isoformas. Las muestras fueron tomadas de sangre y células de granulosa ovárica de mujeres.

Contexto:
El gen ***FMR1*** (ENSG00000102081), localizado en el cromosoma X, codifica la proteína **FMRP** y está involucrado en cuatro trastornos genéticos. A través del *splicing alternativo* de su ARNm, el transcripto puede generar numerosas isoformas, lo que sugiere que cada una podría tener un rol celular específico. Con todos sus exones, sin ningún intrón, y desde el sitio de inicio de la traducción y hasta el sitio de *stop*, el mensajero de ***FMR1*** mide aproximaadmente 3,8 kpb. Los *primers* fueron diseñados para hibridar algunos pb rio arriba del ATG y algunos pb rio abajo del *stop*, por lo que un amplicón del mensajero completo mide un poco más de 3,9 kpb.

En el laboratorio nos interesa la **Insuficiencia ovárica primaria asociada a la fragilidad del X (FXPOI)** y el patrón de expresión de las isoformas en tejido ovárico. En estudios previos, identificamos varias isoformas durante la foliculogénesis en el ovario de rata, pero debido al diseño experimental no pudimos detectar todas las isoformas potencialmente expresadas. De manera similar, no existen estudios en tejidos humanos que describan todas las isoformas expresadas y sus secuencias completas. Por lo tanto, como objetivo, buscamos optimizar la detección de isoformas mediante secuenciación de lecturas largas con tecnología **Oxford Nanopore Technologies (ONT)** en sangre y células de granulosa de mujeres enroladas en un protocolo de ovodonación (no portadoras de la premutación ni la mutación de ***FMR1***). Este enfoque permite secuenciar los transcriptos completos con alta sensibilidad, lo que nos brinda la posibilidad de identificar isoformas nuevas.

---

## 2. Objetivos del trabajo práctico

Al finalizar el TP, se espera que el/la estudiante sea capaz de:

* Comprender el flujo general de análisis de datos de amplicon-seq con ONT
* Alinear lecturas largas contra un genoma de referencia
* Interpretar métricas básicas de alineamiento
* Inferir isoformas transcriptómicas utilizando el pipeline de **FLAIR**
* Comparar isoformas detectadas entre dos muestras
* Visualizar alineamientos e isoformas en un navegador genómico (IGV)
* Entender por qué las lecturas largas exceden las capacidades de tecnologías de lecturas cortas para este tipo de análisis

---

## 3. Organización del directorio de trabajo

El repositorio del TP está organizado de la siguiente manera:

```text
clase_amplicon_seq/
├── ref/
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh38.115.gtf
│   └── isoforms.fa
├── reads/
│   ├── sample_A.fastq
│   └── sample_B.fastq
├── resultados/
│   ├── alineamiento/
│   ├── correct/
│   ├── collapse/
│   └── quantify/
└── scripts/
    └── subsample.sh
```

### Contenido de los directorios

* **ref/**: archivos de referencia (genoma y anotaciones)
* **reads/**: lecturas ONT de cada muestra
* **resultados/**: salidas de cada etapa del pipeline
* **scripts/**: scripts auxiliares (opcional, con fines didácticos)

---

## 4. Datos de referencia

### Genoma

Se utiliza el genoma humano **GRCh38 (Ensembl release 115)**:

* Archivo FASTA: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`

CAMBIAR SEGUN LO QUE DIGA NATY: > Para el TP se trabaja con una versión reducida del genoma, conteniendo únicamente la región correspondiente al gen **FMR1**, con el objetivo de reducir el tamaño de los archivos y el tiempo de cómputo.

### Anotaciones

* Archivo GTF: `Homo_sapiens.GRCh38.115.gtf`

Este archivo contiene las anotaciones génicas y transcriptómicas utilizadas por FLAIR para corregir y colapsar isoformas.

---

DISCLAIMER:

## Control de calidad (QC) de lecturas ONT (paso previo al análisis)

Antes de iniciar cualquier análisis de datos de secuenciación, es fundamental realizar un **control de calidad (Quality Control, QC)** de las lecturas crudas.

En el caso de datos de **Oxford Nanopore Technologies (ONT)**, el QC suele enfocarse principalmente en:

- **Distribución de longitudes de las lecturas**
- **Calidad promedio por lectura (Q-score)**
- Identificación de lecturas truncadas o artefactos

Herramientas comúnmente utilizadas para este paso incluyen, por ejemplo, `NanoPlot`, `NanoQC` o herramientas similares.

### Ejemplos de distribuciones de longitudes

A continuación se muestran ejemplos ilustrativos de histogramas de longitudes de lecturas ONT:

#### Ejemplo 1: muchas lecturas cortas
Este patrón suele indicar:
- fragmentación del cDNA
- problemas en la preparación de la biblioteca
- reads incompletas


![Distribución con lecturas cortas](images/qc_short_reads.png)


#### Ejemplo 2: lecturas del tamaño esperado
Este patrón es el esperado para un experimento de amplicon-seq bien controlado, donde la mayoría de las lecturas tienen una longitud cercana al tamaño del amplicón.


![Distribución del tamaño esperado](images/qc_long_reads.png)


### Filtrado por calidad y longitud

En un pipeline completo, luego del QC se puede aplicar un **filtrado de lecturas**, por ejemplo:

- eliminar lecturas por debajo de un Q-score mínimo
- eliminar lecturas demasiado cortas o demasiado largas respecto al tamaño esperado del amplicón

Este paso permite mejorar la calidad del alineamiento y la detección de isoformas.

### Nota para este trabajo práctico

En este TP **no se realizará el paso de QC ni filtrado**, ya que estos conceptos y herramientas fueron abordados previamente en el TP de **cDNA-seq**.

Para este ejercicio se trabajará directamente con lecturas ya seleccionadas, con el objetivo de focalizarse en:

- alineamiento contra el genoma
- detección de isoformas
- comparación entre muestras



## 5. Alineamiento de lecturas contra el genoma

### ¿Por qué alinear contra el genoma y no contra el transcriptoma?

* Permite detectar isoformas no anotadas
* Evita sesgos introducidos por anotaciones incompletas
* Es un paso necesario para pipelines de detección de isoformas como FLAIR

### Indexado del genoma

El indexado del genoma con **minimap2** consiste en generar una estructura de datos que permite acelerar el alineamiento.

```bash
minimap2 -d Homo_sapiens.GRCh38.dna.primary_assembly.mmi Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

> Para el TP, el genoma ya se provee indexado para ahorrar tiempo.

### Alineamiento de lecturas

Las lecturas ONT se alinean utilizando parámetros *splice-aware*:

```bash
minimap2 -ax splice Homo_sapiens.GRCh38.dna.primary_assembly.mmi reads/sample_A.fastq --splice-flank yes --junc-bonus 10 -o resultados/alineamiento/sample_A.sam
```

El mismo procedimiento se repite para `sample_B.fastq`.

---

## 6. Evaluación del alineamiento

Se utilizan herramientas de **samtools** para evaluar la calidad del alineamiento:

```bash
samtools flagstat resultados/alineamiento/sample_A.sam
```

Este comando resume:

* Número total de lecturas
* Porcentaje de lecturas alineadas
* Lecturas correctamente pareadas (si aplica)

> Discutir en clase qué significa cada métrica y qué se espera en datos de amplicon-seq.

---

## 7. Conversión y procesamiento de formatos

### SAM → BAM ordenado (*sorted*)

```bash
samtools view -bS sample_A.sam | samtools sort -o sample_A.sorted.bam
```

### Indexado del BAM

```bash
samtools index sample_A.sorted.bam
```

### Conversión a BED12

```bash
bamToBed -bed12 -i sample_A.sorted.bam > sample_A.bed
```

> El formato BED12 es requerido por FLAIR y representa explícitamente la estructura exon–intrón de cada lectura.

---

## 8. Detección de isoformas con FLAIR

### FLAIR correct

Corrige los alineamientos utilizando anotaciones conocidas:

```bash
flair correct --query sample_A.bed --genome ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf ref/Homo_sapiens.GRCh38.115.gtf --output resultados/correct/sample_A
```

### FLAIR collapse

Agrupa lecturas corregidas en isoformas únicas:

```bash
flair collapse -g ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf ref/Homo_sapiens.GRCh38.115.gtf -q flair.all_corrected.bed -r reads/sample_A.fastq --output resultados/collapse/sample_A --stringent --check_splice --generate_map
```

### FLAIR quantify

Cuantifica la abundancia de isoformas:

```bash
flair quantify -r reads-manifest.tsv -i flair.collapsed.isoforms.fa --isoform_bed flair.collapsed.isoforms.bed
```

---

## 9. Visualización en IGV

Los archivos `*.sorted.bam`, `*.bai`, `*.bed` y `*.gtf` pueden cargarse en **IGV** para:

* Visualizar alineamientos individuales
* Observar estructuras de exones
* Comparar isoformas entre Sample A y Sample B

Se recomienda navegar a la región del gen ***FMR1*** y discutir diferencias observadas entre las muestras.

---

## 10. Discusión y preguntas

1. ¿Qué isoformas están presentes en ambas muestras?
2. ¿Se detectan isoformas exclusivas de alguna muestra?
3. ¿Qué ventajas ofrece ONT frente a lecturas cortas para este análisis?
4. ¿Qué limitaciones tiene el enfoque de amplicon-seq?

---

## 11. Conclusión

Este TP muestra cómo el uso de lecturas largas permite una caracterización más precisa del transcriptoma, destacando el potencial de ONT para estudios de isoformas que no pueden resolverse adecuadamente con tecnologías de lecturas cortas.
