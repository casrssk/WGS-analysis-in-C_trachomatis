##################################
###fastqc 1
###################################
rule fastqc1:
	input:r1="results/ct_fastq/{sample}_1.fastq.gz", r2="results/ct_fastq/{sample}_2.fastq.gz"
	output:r1="results/fastqc1/{sample}_1_fastqc.html",r2="results/fastqc1/{sample}_2_fastqc.html"
	threads: 30
	log:"logs/{sample}_fastqc1.log"
	conda:wdir + "envs/preprocess.yml"
	shell:
		"""
		fastqc --outdir results/fastqc1/ --extract  -f fastq {input.r1} {input.r2} --thread 16                                                                                                                                                                                                                                                                               
		"""

#################################
#BBnorm
#################################
rule bbnorm:
	input:r1="results/ct_fastq/{sample}_1.fastq.gz", r2="results/ct_fastq/{sample}_2.fastq.gz"
	output:read1 = "results/normalized/{sample}.norm_1.fastq", read2 = "results/normalized/{sample}.norm_2.fastq"
	params:hist_in = "results/normalized/{sample}_input_hist.txt",
	       hist_out = "results/normalized/{sample}_output_hist.txt",
	       out_toss = "results/normalized/{sample}.toss.fastq",
	       passes = "2",
	       options = "ecc=t fixspikes=t prefilter=t"
	threads: 16
	log:"logs/{sample}_bbnorm.log"
	conda: wdir + "envs/preprocess.yml"
	shell:
		"bbnorm.sh in={input.r1} in2={input.r2} hist={params.hist_in} histout={params.hist_out} out={output.read1} out2={output.read2} outt={params.out_toss} passes={params.passes} threads={threads} {params.options}"

##################################
###fastqc 2
###################################
rule fastqc2:
	input:r1="results/normalized/{sample}.norm_1.fastq.gz", r2="results/normalized/{sample}.norm_2.fastq.gz"
	output:r1="results/fastqc2/{sample}.norm_1_fastqc.html",r2="results/fastqc2/{sample}.norm_2_fastqc.html"
	threads: 30
	log:"logs/{sample}_fastqc2.log"
	conda:wdir + "envs/preprocess.yml"
	shell:
		"""
		fastqc --outdir results/fastqc2/ --extract  -f fastq {input.r1} {input.r2} --thread 16
		"""


