############################################
####### Index bam file-samtool
############################################
rule samtools_index:
	input:"results/mapping/addrg/{sample}.addrg.bam"
	output:"results/mapping/addrg/{sample}.addrg.bam.bai"
	priority: 50
	log:"logs/{sample}_ct_dedup_bai.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		samtools index {input} > {output}
		"""

###########################################
############  variant calling
##############################################
rule variant_calling:
	input:bam="results/mapping/addrg/{sample}.addrg.bam", ref="reference/ct_genome_consensus.fasta"
	output:vcf_file="results/var_call/{sample}.vcf"
	log:"logs/{sample}_vcf.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf_file} -ploidy 1 --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
		"""



#########################################
##### normalise and compress the files 
##############################################
rule create_bcfnorm:
	input:vcf_file="results/var_call/{sample}.vcf",ref="reference/ct_genome_consensus.fasta"
	output:bcf_file="results/var_call/{sample}.norm.bcf"
	log:"logs/{sample}_bcf_norm.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		bcftools norm -m-any {input.vcf_file} | bcftools norm -Ob --check-ref w -f {input.ref} >  {output.bcf_file}
		"""

#########################################
#####bcf index- bcftools 
###############################################
rule bcf_index:
	input:bcf_file="results/var_call/{sample}.norm.bcf"
	output:bcf_index="results/var_call/{sample}.norm.bcf.csi"
	priority: 50
	log:"logs/{sample}_bcf_norm_csi.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		bcftools index -f {input.bcf_file} > {output.bcf_index}
		"""




########################################
###########  merge bcf -bcftools
########################################
rule merge_vcfs:
	input:expand("results/var_call/{sample}.norm.bcf",sample=SAMPLES)
	output:"results/var_call/variants.norm.merged.vcf"
	log:"logs/merged_vcf.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		bcftools merge -Ov -m none {input} > {output} 
		"""

#########################################
############  select snps- gatk
##############################################
rule select_snps:
	input:merged="results/var_call/variants.norm.merged.vcf",ref="reference/ct_genome_consensus.fasta"
	output:snps_vcf="results/var_call/snps/snps.vcf"
	log:"logs/snps_vcf.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk SelectVariants -R {input.ref} -V {input.merged} --select-type-to-include SNP -O {output.snps_vcf} 
		"""

#########################################
#############  select indels- gatk
###############################################
rule select_indel:
	input:indel="results/var_call/variants.norm.merged.vcf",ref="reference/ct_genome_consensus.fasta"
	output:indel_vcf="results/var_call/indels/indel.vcf"
	log:"logs/indel_vcf.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		gatk SelectVariants -R {input.ref} -V {input.indel} --select-type-to-include INDEL -O {output.indel_vcf} 
		"""



#########################################
#############  filter snps- gatk
###############################################
rule filter_snps:
	input:filter_snps="results/var_call/snps/snps.vcf",ref="reference/ct_genome_consensus.fasta"
	output:filtered_vcf="results/var_call/filtered/filtered.snps.vcf"
	log:"logs/filtered_vcf.log"
	params:""
	conda: wdir + "envs/environment.yml"
	shell:
		"""
		cat {input.filter_snps} |java -jar 4.3u/snpEff/SnpSift.jar filter " ( QUAL >= 30 ) && (DP >= 10)" > {output.filtered_vcf}
		"""


################################################
############## filter indels- gatk
################################################
rule filter_indel:
	input:filter_indel="results/var_call/indels/indel.vcf",ref="reference/ct_genome_consensus.fasta"
	output:filtered_indel="results/var_call/filtered/filtered.indels.vcf"
	log:"logs/filtered_indels.log"
	params:""
	conda: wdir + "envs/environment.yml"	
	shell:
		"""
		cat {input.filter_indel} |java -jar 4.3u/snpEff/SnpSift.jar filter " ( QUAL >= 30 ) && (DP >= 2)" > {output.filtered_indel}
		"""




