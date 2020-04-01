#######################################
##### edit vcf file
########################################
rule vcf_edit:
	input:"results/var_call/filtered/filtered.snps.vcf"
	output:"results/var_call/filtered/edit.snps.vcf"
	log:"logs/edit_snps.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		cat {input} | sed -e "s/^NC_000117.1/Chromosome/g" -e "s/^ct_genome/Chromosome/g"  > {output}
		"""

#########################################
###########  annotate snps
#############################################
rule snp_eff:
	input:"results/var_call/filtered/edit.snps.vcf"
	output:calls="results/snpeff/snps.annotate.vcf", stats="results/snpeff/snps.html", csvstats="results/snpeff/snps.csv"
	log:"logs/snps_annotation.log"
	params:reference="Chlamydia_trachomatis_d_uw_3_cx"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		java -Xmx4g -jar /Users/satwantkaur/Desktop/project/samples/4.3u/snpEff/snpEff.jar -c /Users/satwantkaur/Desktop/project/samples/4.3u/snpEff/snpEff.config  {params.reference} -s {output.stats} -csvStats {output.csvstats} {input} > {output.calls}
		"""
#######################################
###### edit vcf file for indels
#########################################
rule edit_indel_vcf:
	input:"results/var_call/filtered/filtered.indels.vcf"
	output:"results/var_call/filtered/edit.indels.vcf"
	log:"logs/edit_indels.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		cat {input} | sed -e "s/^NC_000117.1/Chromosome/g" -e "s/^ct_genome/Chromosome/g"  > {output}
		"""


#########################################
############  annotate indels
##############################################
rule indel_eff:
	input:"results/var_call/filtered/edit.indels.vcf"
	output:calls="results/snpeff/indels.annotate.vcf", stats="results/snpeff/indels.html", csvstats="results/snpeff/indels.csv"
	log:"logs/indels_annotation.log"
	params:reference="Chlamydia_trachomatis_d_uw_3_cx"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		java -Xmx4g -jar /Users/satwantkaur/Desktop/project/samples/4.3u/snpEff/snpEff.jar -c /Users/satwantkaur/Desktop/project/samples/4.3u/snpEff/snpEff.config  {params.reference} -s {output.stats} -csvStats {output.csvstats} {input} > {output.calls}
		"""
															
######################################
##vcf-to-tsv for snps
#####################################
rule snps_vcf_to_tsv:
	input:"results/snpeff/snps.annotate.vcf"
	output:"results/var_call/filtered/snps.table"
	log:"logs/snps_table.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F TYPE -F ANN -GF AD -GF DP -GF GT -O {output}
		"""

######################################
###vcf-to-tsv for indels
######################################
rule indels_vcf_to_tsv:
	input:"results/snpeff/indels.annotate.vcf"
	output:"results/var_call/filtered/indels.table"
	log:"logs/indel_table.log"
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F TYPE -F ANN -GF AD -GF DP -GF GT -O {output}
		"""



#########################################
###### multi qc
#########################################
rule multi_qc:
	input:"results/snpeff/indels.csv", 
              "results/snpeff/snps.csv",
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

