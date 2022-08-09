#################################################
# Michael Wang					#
# Meta-scRNA-seq pipeline			#
# https://github.com/fw262/Meta-scRNA-seq	#
#################################################
import pdb
configfile:'config.yaml'
########################################################################################################
# Variables and references
########################################################################################################
MISMATCH = config['GLOBAL']['allowed_mismatch']
DATADIR = config['DATADIR']
STAREXEC = config['STAREXEC']
KRAKEN = config['KRAKEN']
CORES = config['CORES']

rule all:
############## original default call here
	input: expand('{PIPELINE_MAJOR}/{sample}_solo/plots', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])

rule STARsolo:
        input:  read1 = config['DATADIR']+'/{sample}_R1_001.fastq.gz',
                read2 = config['DATADIR']+'/{sample}_R2_001.fastq.gz',
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
		{KRAKEN} --use-names --db {params.krakenDB} <(samtools view -b -f 4 {input} | samtools fasta) | gzip > {output}
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
		kraken='{path}/{sample}_solo/Kraken_barcodes.txt.gz',
	output:
		directory('{path}/{sample}_solo/plots')
	threads: 1
	shell:
		"""
		Rscript performRNAAnalysis.R -m {input.expMat} -k {input.kraken}
		"""
