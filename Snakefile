#################################################
# Roozbeh Abedini Nassab			#
# Drop-Seq Snakemake pipeline			#
# Adapted from the pipeline by Hoohm:		#
# https://github.com/Hoohm/dropseqpipe		#
# And Dr. De Vlaminck's work on it.		#
# last edited July, 4, 2017			#
#################################################
import pdb
########################################################################################################
#/home/fw262/miniconda3/share/dropseq_tools-1.13-0/home/fw262/miniconda3/share/dropseq_tools-1.13-0/home/fw262/miniconda3/share/dropseq_tools-1.13-0# Configfile
########################################################################################################
configfile:'config.yaml'
########################################################################################################
# Variables and references
########################################################################################################
SAMPLEWDIR = config['SAMPWDIR']
MISMATCH = config['GLOBAL']['allowed_mismatch']
DATADIR = config['DATADIR']
GENOMEREF = config['GENOMEREF']
TMPDIR = config['TMPDIR']
gtffile=config['REFGTF']
MERGEBP=str(config['MERGEBP'])
THRESH=str(config['THRESH'])
#FULLSCRIPTPATH=str(config['FULLSCRIPTPATH'])
########################################################################################################
# Executables
########################################################################################################
TMPDIR = config['TMPDIR']
PICARD = config['PICARD']
DROPSEQ = config['DROPSEQ']
STAREXEC = config['STAREXEC']
CORES = config['CORES']
gtfToGenePred=config['GTFTOGENEPRED']
TRINITY_HOME=config["TRINITY_HOME"]
MINIMAP=config["MINIMAP"]
TADPOLE=config["TADPOLE"]

rule all:
############## original default call here
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_kraken.out', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_minimap_mapped.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_tadpole.fa', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	input: expand('{PIPELINE_MAJOR}/{sample}_solo/plots', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}_solo/', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_STARsolo', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: 'kraken_combined.out'

rule STARsolo:
        input:  read1 = config['DATADIR']+'/{sample}_S1_L001_R1_001.fastq.gz',
                read2 = config['DATADIR']+'/{sample}_S1_L001_R2_001.fastq.gz',
                index = config['STAR_IND']
        params: prefix = 'STAR_ind',
                mismatch = MISMATCH,
                BC_start = config['GLOBAL']['BC_range']['first'],
                BC_length = config['GLOBAL']['BC_range']['length'],
                UMI_start = config['GLOBAL']['UMI_range']['first'],
                UMI_length = config['GLOBAL']['UMI_range']['length'],
		read_length = config['GLOBAL']['total_length'],
		tempDir = '{path}/{sample}_STARsolo_tmp',
		outputPref = '{path}/{sample}_solo'
        output: out1='{path}/{sample}_solo/Aligned.sortedByCoord.out.bam'
        threads: CORES
        shell:
              	"""
		rm -rf {params.tempDir}
		\
                {STAREXEC} \
                        --soloType Droplet \
						--soloCBwhitelist None \
						--outSAMtype BAM SortedByCoordinate \
                        --genomeDir {input.index} \
                        --runThreadN {CORES} \
                        --outFilterMismatchNmax={params.mismatch} \
                        --readFilesIn {input.read2} {input.read1} \
                        --readFilesCommand zcat \
                        --genomeLoad NoSharedMemory \
                        --soloCBstart {params.BC_start} \
                        --soloCBlen {params.BC_length} \
                        --soloUMIstart {params.UMI_start} \
                        --soloUMIlen {params.UMI_length} \
						--soloBarcodeReadLength 1 \
						--outSAMattributes NH HI AS nM CB UB \
						--outSAMunmapped Within KeepPairs \
						--outFileNamePrefix {params.outputPref}/
                """

#####################################
# After running standard cell ranger RNA, extract unmapped reads from possorted_bam_unmapped
#####################################
# extract unmapped reads and run Kraken on unmapped reads
rule getUnmappedReads:
	input:
		'{path}/{sample}_solo/Aligned.sortedByCoord.out.bam'
	params:
		krakenDB=config['KRAKEN_DB']
	output:
		temp('{path}/{sample}_solo/possorted_bam_kraken2.out.gz')
	threads: CORES
	shell:
		"""
		kraken2 --use-names --db {params.krakenDB} <(samtools view -b -f 4 {input} | samtools fasta) | gzip > {output}
		"""
# extract cell barcode and UMI of unmapped reads
rule unmappedBarcodes:
	input:
		bam='{path}/{sample}_solo/Aligned.sortedByCoord.out.bam',
		kraken='{path}/{sample}_solo/possorted_bam_kraken2.out.gz'
	output:
		temp('{path}/{sample}_solo/Kraken_barcodes.txt.gz')
	threads: CORES
	shell:
		"""
		paste <(zcat {input.kraken}) <(samtools view -f 4 {input.bam} | grep -o -P '(?<=CB:Z:).*(?=UB:Z:)') <(samtools view -f 4 {input.bam} | grep -o -P '(?<=UB:Z:).*') | gzip > {output}
		"""

# remove duplicate UMIs in column 8
rule removeUMIDuplicates:
	input:
		'{path}/{sample}_solo/Kraken_barcodes.txt.gz'
	output:
		'{path}/{sample}_solo/Kraken_barcodes2.txt.gz'
	threads: CORES
	shell:
		"""
		zcat {input} | awk '!seen[$8]++' | gzip > {output}
		"""
rule gzipCountMatrices:
	input: '{path}/{sample}_solo/Aligned.sortedByCoord.out.bam'
	params: '{path}/{sample}_solo/Solo.out/Gene/filtered'
	output: '{path}/{sample}_solo/Solo.out/Gene/filtered/barcodes.tsv.gz'
	threads: 1
	shell:
		"""
		gzip {params}/*
		"""

rule generateSeuratPlots:
	input:
		expMat='{path}/{sample}_solo/Solo.out/Gene/filtered/barcodes.tsv.gz',
		kraken='{path}/{sample}_solo/Kraken_barcodes2.txt.gz',
	output:
		directory('{path}/{sample}_solo/plots')
	threads: 1
	shell:
		"""
		Rscript performRNAAnalysis.R -m {input.expMat} -k {input.kraken}
		"""
