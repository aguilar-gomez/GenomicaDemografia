# ==============================================================================
# PIPELINE INTRODUCTORIO DE ANÁLISIS GENÓMICO
# Descarga → QC → Limpieza → Mapeo → BAM → VCF
# ==============================================================================

# ==============================================================================
# 0. INSTALAR O DESCARGAR PROGRAMAS NECESARIOS
# ==============================================================================

# Puedes instalar estos programas en Linux o macOS usando conda, apt, o descargando desde sus sitios oficiales.

# BWA (para mapear lecturas)
# Página oficial: http://bio-bwa.sourceforge.net/
# Repositorio GitHub: https://github.com/lh3/bwa
# Instalación (ejemplo con conda): 
conda install -c bioconda bwa

# SAMtools (para manipular archivos SAM/BAM)
# Página oficial: http://www.htslib.org/
# Repositorio GitHub: https://github.com/samtools/samtools
# Instalación: 
conda install -c bioconda samtools

# BCFtools (para llamar variantes y manipular archivos VCF)
# Página oficial: http://www.htslib.org/
# Repositorio GitHub: https://github.com/samtools/bcftools
# Instalación: 
conda install -c bioconda bcftools

# FastQC (para control de calidad de lecturas)
# Página oficial: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Descarga directa (ZIP): https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
# Instalación con conda: 
conda install -c bioconda fastqc

# Trimmomatic (para limpiar lecturas)
# Página oficial: http://www.usadellab.org/cms/?page=trimmomatic
# Descarga directa (JAR): http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
# Instalación con conda: 
conda install -c bioconda trimmomatic

# ------------------------------------------------------------------------------
# 1. DESCARGAR ARCHIVOS DE LECTURAS
# ------------------------------------------------------------------------------
# Ejemplo con datos públicos de NCBI SRA
# Aquí usamos fasterq-dump (parte de SRA Toolkit), o wget si tienes URLs directas

# module load sra-tools
fasterq-dump SRRXXXXXXX -O raw_data

# O usando wget:
wget -O raw_data/sample_R1.fastq.gz "URL_DEL_ARCHIVO_R1"
wget -O raw_data/sample_R2.fastq.gz "URL_DEL_ARCHIVO_R2"

# ------------------------------------------------------------------------------
# 2. CALCULAR FASTQC (EVALUACIÓN DE CALIDAD DE LECTURAS)
# ------------------------------------------------------------------------------
# Genera reportes HTML y ZIP de calidad de secuenciación
fastqc raw_data/*.fastq.gz -o qc_results

# ------------------------------------------------------------------------------
# 3. LIMPIEZA DE LECTURAS (TRIMMING)
# ------------------------------------------------------------------------------
# Ejemplo con Trimmomatic (para reads pareados)
# Elimina adaptadores y bases de baja calidad
trimmomatic PE \
  raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz \
  trimmed/sample_R1_paired.fq.gz trimmed/sample_R1_unpaired.fq.gz \
  trimmed/sample_R2_paired.fq.gz trimmed/sample_R2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# ------------------------------------------------------------------------------
# 4. MAPEAR LECTURAS AL GENOMA DE REFERENCIA (BWA MEM)
# ------------------------------------------------------------------------------
# Primero, indexar el genoma de referencia (solo se hace una vez)
bwa index reference_genome.fasta

# Luego, realizar el alineamiento (produce archivo SAM)
bwa mem -t 4 reference_genome.fasta \
  trimmed/sample_R1_paired.fq.gz trimmed/sample_R2_paired.fq.gz > mapped/sample.sam

# ------------------------------------------------------------------------------
# 5. CONVERTIR SAM A BAM Y ORDENARLO
# ------------------------------------------------------------------------------
samtools view -bS mapped/sample.sam > bam/sample.bam
samtools sort bam/sample.bam -o bam/sample.sorted.bam
samtools index bam/sample.sorted.bam

# ------------------------------------------------------------------------------
# 6. LLAMAR VARIANTES (VCF)
# ------------------------------------------------------------------------------
# Primero, crear archivo de referencia indexado
samtools faidx reference_genome.fasta

# Luego generar el archivo VCF con bcftools
bcftools mpileup -Ou -f reference_genome.fasta bam/sample.sorted.bam | \
bcftools call -mv -Oz -o vcf/sample.vcf.gz

# Indexar el VCF
bcftools index vcf/sample.vcf.gz
