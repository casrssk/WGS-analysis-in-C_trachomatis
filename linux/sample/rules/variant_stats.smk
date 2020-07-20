#######################################
##### edit vcf file
########################################
rule vcf_edit:
	input:"results/var_call/filtered/{sample}.filtered.snps.vcf"
	output:"results/var_call/filtered/{sample}.edit.snps.vcf"
	log:"logs/{sample}_edit_snps.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		cat {input} | sed -e "s/^NC_000117.1/Chromosome/g" -e "s/^ct_genome/Chromosome/g"  > {output}
		"""

#########################################
###########  annotate snps
#############################################
rule snp_eff:
	input:"results/var_call/filtered/{sample}.edit.snps.vcf"
	output:calls="results/snpeff/{sample}.snps.annotate.vcf", stats="results/snpeff/{sample}.snps.html", csvstats="results/snpeff/{sample}.snps.csv"
	log:"logs/{sample}_snps_annotation.log"
	params:reference="Chlamydia_trachomatis_d_uw_3_cx"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		java -Xmx4g -jar 4.3u/snpEff/snpEff.jar -c 4.3u/snpEff/snpEff.config  {params.reference} -s {output.stats} -csvStats {output.csvstats} {input} > {output.calls}
		"""
#######################################
###### edit vcf file for indels
#########################################
rule edit_indel_vcf:
	input:"results/var_call/filtered/{sample}.filtered.indels.vcf"
	output:"results/var_call/filtered/{sample}.edit.indels.vcf"
	log:"logs/{sample}_edit_indels.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		cat {input} | sed -e "s/^NC_000117.1/Chromosome/g" -e "s/^ct_genome/Chromosome/g"  > {output}
		"""


#########################################
############  annotate indels
##############################################
rule indel_eff:
	input:"results/var_call/filtered/{sample}.edit.indels.vcf"
	output:calls="results/snpeff/{sample}.indels.annotate.vcf", stats="results/snpeff/{sample}.indels.html", csvstats="results/snpeff/{sample}.indels.csv"
	log:"logs/{sample}_indels_annotation.log"
	params:reference="Chlamydia_trachomatis_d_uw_3_cx"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		java -Xmx4g -jar 4.3u/snpEff/snpEff.jar -c 4.3u/snpEff/snpEff.config  {params.reference} -s {output.stats} -csvStats {output.csvstats} {input} > {output.calls}
		"""
															
######################################
##vcf-to-tsv for snps
#####################################
rule snps_vcf_to_tsv:
	input:"results/snpeff/{sample}.snps.annotate.vcf"
	output:"results/var_call/filtered/{sample}.snps.table"
	log:"logs/{sample}_snps_table.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk VariantsToTable -V {input} -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -GF AD -GF DP -GF GT -O {output}
		"""

######################################
###vcf-to-tsv for indels
######################################
rule indels_vcf_to_tsv:
	input:"results/snpeff/{sample}.indels.annotate.vcf"
	output:"results/var_call/filtered/{sample}.indels.table"
	log:"logs/{sample}_indel_table.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F TYPE -F ANN -GF AD -GF DP -GF GT -O {output}
		"""



#########################################
###### multi qc
#########################################
rule multi_qc:
	input:expand("results/snpeff/{sample}.indels.csv", sample=SAMPLES),
              expand("results/snpeff/{sample}.snps.csv", sample=SAMPLES),
	      expand("results/mapping/dedup/{sample}_stats/qualimapReport.html", sample=SAMPLES),
              expand("results/mapping/dedup/{sample}.metrics.txt",sample=SAMPLES),
	      expand("results/fastqc2/{sample}.norm_1_fastqc.html", sample=SAMPLES) if config["normalization"]==1 else expand("results/fastqc1/{sample}_1_fastqc.html", sample=SAMPLES), 
	      expand("results/fastqc1/{sample}_1_fastqc.html", sample=SAMPLES),
	output:"results/multiqc/multiqc_report.html"
	params: [config["multiqc_params"]],
	log:"logs/multi_qc"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		multiqc results --outdir results/multiqc 
		"""

