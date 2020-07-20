#################################
#mapping with bowtie -chlamydia
#################################
rule  mapping_ct:
	input:r1="results/normalized/{sample}.norm_1.fastq.gz" if config["normalization"]==1 else "results/ct_fastq/{sample}_1.fastq.gz",
		r2="results/normalized/{sample}.norm_2.fastq.gz" if config["normalization"]==1 else "results/ct_fastq/{sample}_2.fastq.gz"
	output:"results/mapping/{sample}.bam"
	log: "logs/{sample}_ct_bowtie2.log"
	params: 
		index="reference/ct_pan"
	threads:8
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} 2> {log} | \
		samtools view  -bS  > {output}
		"""

#################################
###getting only mapped reads Samtools
###################################
rule map_reads:
	input:"results/mapping/{sample}.bam"
	output:"results/mapping/{sample}.mapped.bam"
	log:"logs/{sample}_mapped.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		samtools view -b -F 4 {input} > {output} 
		"""


#################################
###sorting with Samtools
###################################
rule sorting_reads:
	input:"results/mapping/{sample}.mapped.bam"
	output:"results/mapping/sort/{sample}.sorted.bam"
	log:"logs/{sample}_sorted.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		picard  SortSam I={input} O={output} SORT_ORDER=coordinate
		"""
              

###################################
### mark_duplicates           
###################################
rule mark_duplicates:
	input:"results/mapping/sort/{sample}.sorted.bam"
	output:bam="results/mapping/dedup/{sample}.dedup.bam", metrics="results/mapping/dedup/{sample}.metrics.txt"
	log:"logs/{sample}_ct_dedup.log"
	params:config["Duplicates"],
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		picard MarkDuplicates {params} I={input} O={output.bam} METRICS_FILE={output.metrics} AS=TRUE
		"""

############################################
########### validate_bam
#################################################
rule validate_bam:
	input:"results/mapping/dedup/{sample}.dedup.bam"
	output:"results/mapping/dedup/{sample}_stats/qualimapReport.html"
	log:"logs/{sample}_bam_stats.log"
	params:outdir="results/mapping/dedup/{sample}_stats"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		qualimap bamqc -bam {input} --outdir {params.outdir} 
		"""

########################################
#### Add or replace read groups-picard             
#########################################
rule addorreplacgroups:
	input:"results/mapping/dedup/{sample}.dedup.bam"
	output:"results/mapping/addrg/{sample}.addrg.bam"
	log:"logs/{sample}_ct_dedup_addrg.log"
	params:config["addorreplacgroups"],
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		picard AddOrReplaceReadGroups I={input} O={output} {params}
		"""

		



