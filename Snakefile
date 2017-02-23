
configfile: "config.yml"
workdir: config["workdir"]

# Decompress fastq.gz 
rule unzip_fastq:
	input:
		"{filename}.fastq.gz"
	output:
		"{filename}.fastq"
	shell:
		"gzip -d {input}"

# Aligning fastq to bampwd
rule alignement: 
	input: 
		forward="{sample}_1.fastq",
		reverse="{sample}_2.fastq"

	output:
		"{sample}.sam"
	threads:128 
	shell:
		"bwa mem -t {threads} {config[ref]} {input.forward} {input.reverse} > {output}"

# Aligning fastq to bam
# rule alignement: 
# 	input: 
# 		forward="{sample}_1.fastq",
# 		reverse="{sample}_2.fastq"

# 	output:
# 		"{sample}.sam"
# 	threads:128 
# 	shell:
# 		"bowtie2 -x /PUBLIC_DATA/ReferenceGenomes/IonTorrenthg19/hg19  -1 {input.forward} -2 {input.reverse} -S {output}"


# Rule sam to bam 
rule sam_to_bam:
	input:
		"{sample}.sam"
	output:
		bam   = "{sample}.bam",
		index = "{sample}.bam.bai" 
	shell:
		"samtools view  -Sb {input} | samtools sort - > {output.bam} && samtools index {output.bam} {output.index}"