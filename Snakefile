
# HugoDims2
# Copyright (C) 2017  Sacha Schutz

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
#I would prefere to use : IlluminaBasecallsToSam
rule alignement: 
	input: 
		forward="{sample}_1.fastq",
		reverse="{sample}_2.fastq"

	output:
		"{sample}.sam"
	threads:128 
	shell:
		"bwa mem -t {threads} {config[ref]} {input.forward} {input.reverse} > {output}"

rule sort:
	input:
		"{filename}.sam"
	output:
		"{filename}.sorted.sam"
	shell:
		"picard-tools SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate"


rule mark_duplicate:
	input:
		"{sample}.sorted.sam"
	output:
		"{sample}.markDuplicate.sam",
		"{sample}.markDuplicate.metrics",
	shell:
		"picard-tools MarkDuplicatesWithMateCigar INPUT={input} OUTPUT={output[0]} M={output[1]}"


#Rule sam to bam 
rule sam_to_bam:
	input:
		"{filename}.sam"
	output:
		bam   = "{filename}.bam",
	shell:
		"samtools view  -Sbh {input} | samtools sort - > {output.bam}"

rule index_bam:
	input:
		"{filename}.bam"
	output:
		"{filename}.bam.bai"
	shell:
		"samtools index {input} {output}"


rule add_read_group:
	input:
		"{sample}.bam",
		"{sample}.bam.bai"
	output:
		"{sample}.group.bam"
	shell:
		"picard-tools AddOrReplaceReadGroups I={input[0]} O={output} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"




#Haplotype caller 

rule haplotype_caller:
	input:
		"{sample}.group.bam",
		"{sample}.group.bam.bai"
	output:
		"{sample}.g.vcf"
	shell:
		"gatk -T HaplotypeCaller -R {config[ref]} -I {input[0]} -o {output} -ERC GVCF"



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


