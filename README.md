# GenomicaDemografia
Taller de Análisis Genómicos de Diversidad Genética e Historia Demográfica de Poblaciones de Mamíferos Acuáticos

## 1. Formatos de Archivos Genómicos

En esta sección aprenderemos los conceptos básicos necesarios para iniciar un análisis genómico.
Nos centraremos en tres objetivos principales:

- Conocer las tecnologías de secuenciación más utilizadas en estudios de diversidad genética e historia demográfica, y cómo sus características influyen en los tipos de datos que obtenemos.

- Familiarizarnos con los formatos de archivo más comunes en bioinformática (FASTQ, FASTA, SAM/BAM, VCF), comprendiendo qué información contiene cada uno y en qué etapa del análisis se generan o utilizan.

- Seguir un flujo de trabajo introductorio (pipeline) para generar estos archivos paso a paso, desde la descarga de lecturas crudas hasta la obtención de un archivo de variantes (VCF), utilizando herramientas ampliamente empleadas como FastQC, Trimmomatic, BWA, SAMtools y BCFtools.

Esta primera parte sentará las bases para comprender los procesos posteriores del taller, donde analizaremos la calidad de los datos, el alineamiento al genoma de referencia y la detección de variantes.


## 3. Estructura Poblacional: PCA y ADMIXTURE

En esta sección exploraremos cómo inferir y visualizar la estructura genética entre individuos o poblaciones utilizando dos métodos ampliamente empleados en genómica de poblaciones:

PCA (Análisis de Componentes Principales), que resume la variación genética en un espacio de pocas dimensiones y permite identificar agrupamientos, gradientes de variación y posibles individuos atípicos o mezclados.

ADMIXTURE, un algoritmo de clustering basado en modelos que estima la proporción de ancestrías de cada individuo bajo distintos valores de K (número de grupos genéticos), ayudándonos a interpretar procesos como mezcla reciente, diferenciación poblacional y aislamiento histórico.

Además, aprenderemos a preparar los datos correctamente para estos análisis:

- Filtrar variantes para quedarnos con SNPs bialélicos de alta calidad.

- Aplicar filtros opcionales como MAF, thinning, y remoción de sitios con linkage para obtener un conjunto de marcadores independientes.

- Convertir el VCF a formato PLINK y generar los archivos necesarios para cada método.

Esta parte nos permitirá identificar patrones de estructura poblacional en nuestros datos, visualizar relaciones entre muestras y establecer el contexto para análisis más avanzados de demografía e historia evolutiva.
