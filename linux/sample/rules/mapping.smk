##################################################################
### mapping reads to human reference genome usinhg bowtie- human
###################################################################
rule bowtie2:
	input:r1=wdir + "raw/{sample}_1.fastq", r2=wdir + "raw/{sample}_2.fastq"
	output:"bam_files/{sample}.bam"
	log: "logs/{sample}_bowtie2.log"
	params: 
		index="reference/human_genome/GRCh38"
	threads:8
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} 2> {log} | \
		samtools view  -bS  > {output}
		"""

##############################################################
### Filter unmapped reads to get the chlamydia related reads              
##############################################################
rule filter:
	input:"bam_files/{sample}.bam"
	output:"bam_files/{sample}.unmappedreads.bam"
	log:"logs/{sample}.unmappedreads.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		samtools view -b -f 12 -F 256 {input} > {output}
		"""


##############################################
### sorting the unmapped reads with Samtools
###############################################
rule sorting:
	input:"bam_files/{sample}.unmappedreads.bam"
	output:"bam_files/{sample}.sorted.bam"
	log:"logs/{sample}_sorted.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		samtools sort -n {input} > {output}
		"""

############################################################
### converting the unmapped sorted bam files to fastq files
############################################################
rule samtools_bam2fq_interleaved:
	input:"bam_files/{sample}.sorted.bam"
	output:R1="results/ct_fastq/{sample}_1.fastq", R2="results/ct_fastq/{sample}_2.fastq"
	params:" "
	threads: 2
	conda: wdir + "envs/environment.yml"
	shell:
        	"bedtools bamtofastq -i {input} -fq {output.R1} -fq2 {output.R2}"

