import os
import pandas as pd
import tempfile

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: "config.json"

TMPDIR = tempfile.gettempdir()
THREADS = str(config.get("threads", 1))
CHROMOSOMES = ("chr1", "chr2", "chr3", "chr4", "chr5")

PARTS = int(config.get("parts", 10))
manifest = config["manifest"]
samples = config["samples"]
instruments = ["Illumina Genome Analyzer",
               "Illumina Genome Analyzer II",
               "Illumina Genome Analyzer IIx",
               "Illumina HiSeq 2000"]
columns = ("sample_accession", "experiment_accession", "run_accession", "fastq_ftp", "instrument_model")
runs_by_sample_and_experiment = {}
urls_by_run = {}
df = pd.read_table(manifest)

for record in df.ix[(df["sample_accession"].isin(samples)) & (df["instrument_model"].isin(instruments)) & (df["library_layout"] == "PAIRED"), columns].values:
    # Setup runs by sample and experiment.
    if not record[0] in runs_by_sample_and_experiment:
        runs_by_sample_and_experiment[record[0]] = {}

    if not record[1] in runs_by_sample_and_experiment[record[0]]:
        runs_by_sample_and_experiment[record[0]][record[1]] = []

    runs_by_sample_and_experiment[record[0]][record[1]].append(record[2])

    # Setup urls by run.
    if not record[0] in urls_by_run:
        urls_by_run[record[0]] = {}

    if not record[1] in urls_by_run[record[0]]:
        urls_by_run[record[0]][record[1]] = {}

    # urls_by_run: {sample: {experiment: {run1: [run1.1.fastq.gz, run1.2.fastq.gz]}}}
    urls_by_run[record[0]][record[1]][record[2]] = [url.replace("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:") for url in record[3].split(";")]

def _get_bwa_output_for_runs_by_sample(wildcards):
    return ["alignments/%s/%s/%s/alignments.bam" % (wildcards.sample, experiment, run)
            for experiment in runs_by_sample_and_experiment[wildcards.sample]
            for run in runs_by_sample_and_experiment[wildcards.sample][experiment]]

def _get_mrsfast_output_for_runs_by_sample(wildcards):
    outputs = []
    for experiment in runs_by_sample_and_experiment[wildcards.sample]:
        runs = runs_by_sample_and_experiment[wildcards.sample][experiment]
        for run in runs:
            for i in range(len(urls_by_run[wildcards.sample][experiment][run])):
                outputs.append("mrsfast_depths/%s/%s_%s_%s.depth" % (wildcards.sample, experiment, run, i + 1))

    return outputs

#list(runs_by_sample_and_experiment.keys())

rule all:
    input:
        expand("stats_plots/{sample}/gc-depth.png", sample=samples),
        expand("tracks/cnvs.{sample}.bw", sample=samples),
        expand("cnvs_compressed/{sample}.bed.gz", sample=samples),
        expand("merged_chunked_mrsfast_alignments/{sample}.bam", sample=samples),
        expand("smoothed_read_depth_bigWig/{sample}.bw", sample=samples)

#
# Summarize alignments.
#

rule convert_excess_depth_regions_to_excluded_intervals_for_GATK:
    input: "smoothed_depth_gt40.bed"
    output: "excluded_regions_for_GATK.intervals"
    conda: "envs/snake-cnv.yaml"
    shell: """awk '{{ print $1":"$2"-"$3 }}' {input} > {output}"""

rule plot_copy_number_status_for_cnvs_for_all_samples:
    input: "copy_number_in_large_cnvs.high_quality.tab"
    output: "copy_number_in_large_cnvs.high_quality.pdf"
    conda: "envs/snake-cnv.yaml"
    shell: "Rscript plot_copy_number_in_large_cnvs.R {input} {output}"

rule filter_copy_number_status_for_cnvs_by_high_quality_samples:
    input: "high_quality_samples_by_Tgt_copy_number.txt", "copy_number_in_large_cnvs.tab"
    output: "copy_number_in_large_cnvs.high_quality.tab"
    conda: "envs/snake-cnv.yaml"
    shell: "head -n 1 {input[1]} > {output}; grep -wf {input} >> {output}"

rule collect_copy_number_status_for_cnvs:
    input: expand("copy_number_in_large_cnvs/{sample}.tab", sample=list(runs_by_sample_and_experiment.keys()))
    output: "copy_number_in_large_cnvs.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """head -n 1 {input[0]} | awk 'OFS="\\t" {{ print $0,"sample" }}' > {output}; for file in {input}; do sample=`basename $file | sed 's/.tab//'`; sed 1d $file | awk -v sample="$sample" 'OFS="\\t" {{ print $0,sample }}'; done >> {output}"""

rule filter_copy_number_by_region_by_high_quality_samples:
    input: "high_quality_samples_by_Tgt_copy_number.txt", "copy_number_by_region.tab"
    output: "copy_number_by_region.high_quality.tab"
    conda: "envs/snake-cnv.yaml"
    shell: "head -n 1 {input[1]} > {output}; grep -wf {input} >> {output}"

rule get_high_quality_samples:
    input: "copy_number_by_region.tab"
    output: "high_quality_samples_by_Tgt_copy_number.txt"
    conda: "envs/snake-cnv.yaml"
    shell: """awk '($1 == "Tgt" || $1 == "Dp") && ($3 <= 2.5 && $3 >= 1.5)' {input} | cut -f 2 | sort | uniq > {output}"""

rule calculate_copy_number_across_regions:
    input: expand("cnvs_compressed/{sample}.bed.gz", sample=list(runs_by_sample_and_experiment.keys()))
    output: "copy_number_by_region.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """touch cnvs_compressed/*.tbi; cat regions.tab | while read line; do set -- $line; region_name=$1; region=$2; for file in {input}; do sample=`basename $file | sed 's/.bed.gz//'`; tabix $file $region | awk 'OFS="\\t" {{ bases = $3 - $2; print $0,bases,bases*$4 }}' | bedtools groupby -i stdin -g 1 -c 5,6 -o sum,sum | awk -v region=$region_name -v sample=$sample 'OFS="\\t" {{ print region,sample,$3 / $2 }}'; done; done | awk 'OFS="\\t" {{ if (NR == 1) {{ print "region","sample","copy_number" }} print }}' > {output}"""

rule plot_copy_2_status_for_cnvs:
    input: "copy_number_by_diploid_status/{sample}.tab"
    output: "copy_number_by_diploid_status_plots/{sample}.pdf"
    conda: "envs/snake-cnv.yaml"
    shell: "Rscript plot_copy_number_by_status.R {input} {output}"

rule plot_copy_number_status_for_cnvs:
    input: "copy_number_in_large_cnvs/{sample}.tab"
    output: "copy_number_in_large_cnvs_plots/{sample}.pdf"
    conda: "envs/snake-cnv.yaml"
    shell: "Rscript plot_copy_number_in_large_cnvs.R {input} {output}"

rule collect_copy_numbers_by_dels_and_dups:
    input: "copy_number_in_large_duplications/{sample}.tab", "copy_number_in_large_deletions/{sample}.tab", "copy_number_in_copy_2_regions/{sample}.tab"
    output: "copy_number_in_large_cnvs/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """awk 'OFS="\\t" {{ if (NR == 1) {{ print "mean_cn","type" }} print }}' {input} > {output}"""

rule calculate_copy_numbers_in_copy_2_regions:
    input: cnvs="cnvs_bed/{sample}.bed", copy_2="annotations/filtered_copy_2_regions.bed"
    output: "copy_number_in_copy_2_regions/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """bedtools intersect -a {input.copy_2} -b {input.cnvs} -wo | awk 'OFS="\\t" {{ print $1,$2,$3,$7,$8,sprintf("%i", $7*$8) }}' | bedtools groupby -i stdin -g 1,2,3 -c 5,6 -o sum,sum | awk 'OFS="\\t" {{ print $5/$4,"copy_2" }}' > {output}"""

rule get_copy_2_regions_without_deletions_or_duplications:
    input: copy_2="annotations/copy_2_regions.bed", duplications="annotations/duplicated_regions.bed", deletions="annotations/large_deletions.bed"
    output: "annotations/filtered_copy_2_regions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "sort -k 1,1 -k 2,2n {input.duplications} {input.deletions} | bedtools merge -i stdin -d 0 | bedtools intersect -a {input.copy_2} -b stdin -v -sorted > {output}"

rule calculate_copy_numbers_in_large_duplications:
    input: cnvs="cnvs_bed/{sample}.bed", duplications="annotations/duplicated_regions.bed"
    output: "copy_number_in_large_duplications/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """bedtools intersect -a {input.duplications} -b {input.cnvs} -wo | awk 'OFS="\\t" {{ print $1,$2,$3,$7,$8,sprintf("%i", $7*$8) }}' | bedtools groupby -i stdin -g 1,2,3 -c 5,6 -o sum,sum | awk 'OFS="\\t" {{ print $5/$4,"duplication" }}' > {output}"""

rule calculate_copy_numbers_in_large_deletions:
    input: cnvs="cnvs_bed/{sample}.bed", deletions="annotations/large_deletions.bed"
    output: "copy_number_in_large_deletions/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """bedtools intersect -a {input.deletions} -b {input.cnvs} -wo | awk 'OFS="\\t" {{ print $1,$2,$3,$7,$8,sprintf("%i", $7*$8) }}' | bedtools groupby -i stdin -g 1,2,3 -c 5,6 -o sum,sum | awk 'OFS="\\t" {{ print $5/$4,"deletion" }}' > {output}"""

rule annotate_copy_2_status_for_cnvs:
    input: cnvs="cnvs_bed/{sample}.bed", copy_2="annotations/copy_2_regions.bed"
    output: "copy_number_by_diploid_status/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """bedtools intersect -a {input.cnvs} -b {input.copy_2} -wao -sorted | awk 'OFS="\\t" {{ if (NR == 1) {{ print "copy_number","diploid_status" }} if ($5 == ".") {{ type="not_cp2" }} else {{ type="cp2" }} print $4,type }}' > {output}"""

rule create_cnvs_bigwig:
    input: "cnvs_bed/{sample}.bed", "dm6.chrom.sizes"
    output: "tracks/cnvs.{sample}.bw"
    conda: "envs/snake-cnv.yaml"
    shell: "bedGraphToBigWig {input} {output}"

rule compress_and_index_cnvs_bed:
    input: "cnvs_bed/{sample}.bed"
    output: "cnvs_compressed/{sample}.bed.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "bgzip -c {input} > {output}; tabix -f -p bed {output}"

rule prepare_cnvs_bed:
    input: "cnvs/{sample}.copynumber.bed"
    output: "cnvs_bed/{sample}.bed"
    conda: "envs/snake-cnv.yaml"
    shell: """sed '1,2d' {input} | sort -k 1,1 -k 2,2n | cut -f 1-3,5 | awk 'OFS="\\t" {{ copy = sprintf("%.0f", $4); print $1,$2,$3,$4 }}' > {output}"""

rule call_cnvs:
    input: configuration="reference/dm6.cnvr", depths="chunked_mrsfast_depths/{sample}.depth"
    output: "cnvs/{sample}.copynumber.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "mrcanavar --call --xx --multgc -conf {input.configuration} -depth {input.depths} -o cnvs/{wildcards.sample}"

rule merge_depths_by_sample:
    input: configuration="reference/dm6.cnvr", depths=_get_mrsfast_output_for_runs_by_sample
    output: "merged_mrsfast_depths/{sample}.depth"
    conda: "envs/snake-cnv.yaml"
    shell: "mrcanavar --conc -conf {input.configuration} -concdepth {input.depths} -depth {output}"

rule calculate_depths_from_chunked_alignments:
    input: configuration="reference/dm6.cnvr", alignments=expand("chunked_mrsfast_alignments/{{sample}}/{part}.sam.gz", part=range(0, PARTS + 1))
    output: "chunked_mrsfast_depths/{sample}.depth"
    conda: "envs/snake-cnv.yaml"
    shell: "mrcanavar --read -conf {input.configuration} -samdir `dirname {input.alignments}` -depth {output} --gz"

rule convert_depth_bedGraph_to_bigWig:
    input: "smoothed_read_depth_bedGraph/{sample}.bedGraph", "reference/dm6.fasta.fai"
    output: "smoothed_read_depth_bigWig/{sample}.bw"
    conda: "envs/snake-cnv.yaml"
    shell: "bedGraphToBigWig {input} {output}"

rule convert_smoother_output_to_bedGraph:
    input: "smoothed_read_depth/{sample}.bed"
    output: "smoothed_read_depth_bedGraph/{sample}.bedGraph"
    conda: "envs/snake-cnv.yaml"
    shell: """awk 'OFS="\\t" {{ print $1,$2,$2 + 10,$5 }}' {input} > {output}"""

rule get_smoothed_read_depth:
    input: "merged_chunked_mrsfast_depth/{sample}.tab"
    output: "smoothed_read_depth/{sample}.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "smoother -o col3 -w 50 -s 10 -t -f {input} > {output}"

rule get_depth_for_merged_chunked_reads:
    input: "merged_chunked_mrsfast_alignments/{sample}.bam"
    output: "merged_chunked_mrsfast_depth/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: "samtools depth -aa {input} > {output}"

rule merge_chunked_read_bams:
    input: expand("sorted_chunked_mrsfast_alignments/{{sample}}/{part}.bam", part=range(0, PARTS + 1))
    output: "merged_chunked_mrsfast_alignments/{sample}.bam"
    conda: "envs/snake-cnv.yaml"
    shell: "samtools merge {output} {input}"

rule convert_mapped_chunked_reads_from_sam_to_bam:
    input: "chunked_mrsfast_alignments/{sample}/{part}.sam.gz"
    output: temp("sorted_chunked_mrsfast_alignments/{sample}/{part}.bam")
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | samtools sort > {output}"

#
# Align reads with mrsFAST-ULTRA
#

rule map_chunked_reads_with_mrsfast:
    input: reads="merged_alignments/{sample}.bam", index="reference/dm6.masked.fasta.index", reference="reference/dm6.masked.fasta", snps="dgrp2_dm6_dbSNP.index"
    output: "chunked_mrsfast_alignments/{sample}/{part}.sam.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "bam_chunker_cascade -b {input.reads} -p {wildcards.part} -n {PARTS} -u 1 2>> /dev/stderr | mrsfast --search {input.reference} --seq /dev/stdin -o `echo {output} | sed 's/.sam.gz//'` --snp {input.snps} --outcomp -e 2 --crop 36 --threads=1 --mem 6 -n 0  --disable-nohits"

rule index_snps_for_mrsfast:
    input: "dgrp2_dm6_dbSNP.vcf"
    output: "dgrp2_dm6_dbSNP.index"
    conda: "envs/snake-cnv.yaml"
    shell: "snp_indexer {input} {output}"

rule uncompress_dgrp_variants:
    input: "dgrp2_dm6_dbSNP.vcf.gz"
    output: "dgrp2_dm6_dbSNP.vcf"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} > {output}"

rule download_dgrp_variants:
    output: "dgrp2_dm6_dbSNP.vcf.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "wget https://zenodo.org/record/155396/files/{output} -O {output}"

rule index_reference_for_mrsfast:
    input: "reference/dm6.masked.fasta"
    output: "reference/dm6.masked.fasta.index"
    conda: "envs/snake-cnv.yaml"
    shell: "mrsfast --index {input}"

rule prepare_reference_for_mrcanavar:
    input: reference="reference/dm6.masked.fasta", gaps="annotations/gaps.bed"
    output: "reference/dm6.cnvr"
    params: large_window=1000, small_window=500
    conda: "envs/snake-cnv.yaml"
    shell: "mrcanavar --prep -fasta {input.reference} -gaps {input.gaps} -conf {output} -lw_size {params.large_window} -lw_slide {params.small_window} -sw_size {params.small_window} -sw_slide {params.small_window} -cw_size {params.small_window}"

rule hardmask_reference:
    input: reference="reference/dm6.unmasked.fasta", repeats="dm6.repeats_without_RNA.36bp_padded.bed"
    output: "reference/dm6.masked.fasta"
    conda: "envs/snake-cnv.yaml"
    shell: "bedtools maskfasta -fi {input.reference} -bed {input.repeats} -fo {output}"

rule unmask_reference:
    input: "reference/dm6.fasta"
    output: "reference/dm6.unmasked.fasta"
    conda: "envs/snake-cnv.yaml"
    shell: "sed '/^>/!s/\(.*\)/\\U\\1/' {input} > {output}"

rule find_gaps_in_reference:
    input: "reference/dm6.fasta"
    output: "annotations/gaps.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "python scripts/find_fasta_gaps.py {input} | bedtools merge -i stdin -d 0 > {output}"

#
# Align reads with BWA MEM.
#

rule plot_alignment_stats:
    input: "stats/{sample}.tab"
    output: "stats_plots/{sample}/gc-depth.png"
    conda: "envs/snake-cnv.yaml"
    shell: "plot-bamstats -p stats_plots/{wildcards.sample}/ {input}"

rule get_alignment_stats:
    input: "merged_alignments/{sample}.bam"
    output: "stats/{sample}.tab"
    conda: "envs/snake-cnv.yaml"
    shell: "samtools stats {input} > {output}"

rule merge_alignments_by_sample:
    input: _get_bwa_output_for_runs_by_sample
    output: protected("merged_alignments/{sample}.bam")
    conda: "envs/snake-cnv.yaml"
    shell: "samtools merge {output} {input}; samtools index {output}"

rule map_reads:
    input:
        reference="reference/dm6.fasta",
        bwa_index="reference/dm6.fasta.bwt",
        reads_1="reads/{sample}/{experiment}/{run}/1.fastq.gz",
        reads_2="reads/{sample}/{experiment}/{run}/2.fastq.gz"
    output:
        alignments=temp("alignments/{sample}/{experiment}/{run}/alignments.bam"),
        discordant="alignments/{sample}/{experiment}/{run}/discordant.sam",
        split="alignments/{sample}/{experiment}/{run}/split.sam"
    params: threads="2"
    conda: "envs/snake-cnv.yaml"
    shell: "bwa mem -R '@RG\\tID:{wildcards.sample}_{wildcards.run}\\tSM:{wildcards.sample}\\tLB:{wildcards.experiment}\\tPL:illumina\\tPU:{wildcards.run}' -t {params.threads} {input.reference} {input.reads_1} {input.reads_2} | samtools view -h -S -f 1 - | samblaster --ignoreUnmated -r -e -d {output.discordant} -s {output.split} | samtools sort -O bam -T $TMPDIR/jlhudd_{wildcards.run} > {output.alignments}"

#
# Download reads.
#

def _get_url_by_run_and_number(wildcards):
    return urls_by_run[wildcards.sample][wildcards.experiment][wildcards.run][int(wildcards.number) - 1]

rule download_reads_by_run_and_number:
    output: temp("reads/{sample}/{experiment}/{run}/{number}.fastq.gz")
    params:
        url=_get_url_by_run_and_number,
        aspera_key_path=config["aspera_key_path"]
    conda: "envs/snake-cnv.yaml"
    shell: """ascp -QT -l 300m -P33001 -i {params.aspera_key_path} {params.url} {output}"""

#
# Get expected copy number 2 regions
#

rule merge_cnv_duplications_and_segmental_duplications:
    input: "annotations/segmental_duplications.dm6.sorted.bed", "annotations/cnv_duplications.bed"
    output: "annotations/duplicated_regions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "sort -k 1,1 -k 2,2n {input} | cut -f 1-3 | bedtools merge -i stdin -d 0 | awk '$3 - $2 > 1000' > {output}"

rule get_common_large_emerson_duplications:
    input: "../../data/emerson_cnvs.bed"
    output: "annotations/cnv_duplications.bed"
    conda: "envs/snake-cnv.yaml"
    shell: """awk '$4 == "dup" && $3 - $2 > 500 && $5 >= 8' {input} | cut -f 1-3 | bedtools merge -i stdin -d 0 > {output}"""

rule get_copy_2_regions:
    input: bed="annotations/non_copy_2_regions.1kbp_slop.bed", chromosomes="dm6.chrom.sizes"
    output: "annotations/copy_2_regions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "bedtools complement -i {input.bed} -g {input.chromosomes} | sort -k 1,1 -k 2,2n > {output}"

rule add_slop_to_non_copy_2_regions:
    input: bed="annotations/non_copy_2_regions.bed", chromosomes="dm6.chrom.sizes"
    output: "annotations/non_copy_2_regions.1kbp_slop.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "bedtools slop -i {input.bed} -g {input.chromosomes} -b 1000 | bedtools merge -i stdin -d 0 > {output}"

rule merge_cnvs_and_segmental_duplications:
    input: "annotations/segmental_duplications.dm6.sorted.bed", "../../data/emerson_cnvs.bed"
    output: "annotations/non_copy_2_regions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "sort -k 1,1 -k 2,2n {input} | cut -f 1-3 | bedtools merge -i stdin -d 0 > {output}"

#
# Prepare DGRP SVs
#

rule merge_deletions_from_dgrp_and_emerson:
    input: "annotations/dgrp2_large_deletions.bed", "annotations/emerson_common_large_deletions.bed"
    output: "annotations/large_deletions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "cut -f 1-3 {input} | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d 0 > {output}"

rule get_large_dgrp_deletions:
    input: "dgrp2_sv.dm6.bed"
    output: "annotations/dgrp2_large_deletions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: """awk '$4 == "del" && $3 - $2 > 1000' {input} > {output}"""

rule get_large_emerson_deletions:
    input: "../../data/emerson_cnvs.bed"
    output: "annotations/emerson_common_large_deletions.bed"
    conda: "envs/snake-cnv.yaml"
    shell: """awk '$4 == "del" && $3 - $2 > 500 && $5 >= 8' {input} > {output}"""

rule convert_dgrp_all_calls_to_bigbed:
    input: "dgrp2_complete.dm6.bed", "dm6.chrom.sizes"
    output: "tracks/dgrp2_complete.dm6.bb"
    conda: "envs/snake-cnv.yaml"
    shell: "bedToBigBed {input} {output}"

rule liftover_dgrp_all_calls_for_dm3_to_dm6:
    input: "dgrp2_complete.dm3.bed", "dm3ToDm6.over.chain"
    output: "dgrp2_complete.dm6.bed", "dgrp2_complete.unlifted.txt"
    conda: "envs/snake-cnv.yaml"
    shell: "liftOver -minMatch=0.95 {input} {output}; sort -k 1,1 -k 2,2n -o {output[0]} {output[0]}"

rule get_all_dgrp_calls:
    input: "dgrp2.vcf"
    output: "dgrp2_complete.dm3.bed"
    conda: "envs/snake-cnv.yaml"
    shell: """sed '/^#/d' {input} | cut -f 1-5 | awk 'OFS="\\t" {{ ref=length($4); alt=length($5); if (ref - alt > 1) {{ type = "del"; size = ref - alt }} else if (alt - ref > 1) {{ type = "ins"; size = 1 }} else {{ type = "snv"; size = 1 }} print "chr"$1,$2 - 1,$2 - 1 + size,type }}' | sort -k 1,1 -k 2,2n | uniq > {output}"""

rule liftover_dgrp_sv_calls_for_dm3_to_dm6:
    input: "dgrp2_sv.dm3.bed", "dm3ToDm6.over.chain"
    output: "dgrp2_sv.dm6.bed", "dgrp2_sv.unlifted.txt"
    conda: "envs/snake-cnv.yaml"
    shell: "liftOver -minMatch=0.95 {input} {output}"

rule get_dgrp_sv_calls:
    input: "dgrp2.vcf"
    output: "dgrp2_sv.dm3.bed"
    conda: "envs/snake-cnv.yaml"
    shell: """sed '/^#/d;/SNP/d' {input} | cut -f 1-5 | awk 'OFS="\\t" {{ ref=length($4); alt=length($5); if (ref - alt > 50) {{ type = "del"; size = ref - alt }} else if (alt - ref > 50) {{ type = "ins"; size = 1 }} else {{ type = "" }} if (type == "del" || type == "ins") {{ print "chr"$1,$2 - 1,$2 - 1 + size,type }} }}' | sort -k 1,1 -k 2,2n | uniq > {output}"""

rule download_complete_dgrp_calls:
    output: "dgrp2.vcf"
    conda: "envs/snake-cnv.yaml"
    shell: """wget "http://dgrp2.gnets.ncsu.edu/data/website/{output}" -O {output}"""

#
# Prepare segmental duplications.
#

rule build_bigbed_for_duplications:
    input: "annotations/segmental_duplications.dm6.sorted.bed", "dm6.chrom.sizes"
    output: "tracks/segmental_duplications.bb"
    conda: "envs/snake-cnv.yaml"
    shell: "bedToBigBed {input} {output}"

rule merge_duplications_in_dm6:
    input: "annotations/segmental_duplications.dm6.sorted.bed"
    output: "annotations/segmental_duplications.dm6.merged.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "bedtools merge -i {input} -d 0 > {output}"

rule sort_duplications_in_dm6:
    input: "annotations/segmental_duplications.dm6.bed"
    output: "annotations/segmental_duplications.dm6.sorted.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "sort -k 1,1 -k 2,2n {input} > {output}"

rule liftover_segmental_duplications_for_dm3_to_dm6:
    input: "annotations/segmental_duplications.dm3.bed", "dm3ToDm6.over.chain"
    output: "annotations/segmental_duplications.dm6.bed", "annotations/segmental_duplications.unlifted.txt"
    conda: "envs/snake-cnv.yaml"
    shell: "liftOver -minMatch=0.75 {input} {output}"

rule uncompress_liftover_for_dm3_to_dm6:
    input: "dm3ToDm6.over.chain.gz"
    output: "dm3ToDm6.over.chain"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} > {output}"

rule download_liftover_for_dm3_to_dm6:
    output: "dm3ToDm6.over.chain.gz"
    conda: "envs/snake-cnv.yaml"
    shell: """wget "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/liftOver/{output}" -O {output}"""

rule convert_segmental_duplications_for_dm3_to_bed:
    input: "annotations/segmental_duplications.dm3.tab"
    output: "annotations/segmental_duplications.dm3.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "cut -f 1-3 {input} | sort -k 1,1 -k 2,2n > {output}"

rule download_segmental_duplications_for_dm3:
    output: "annotations/segmental_duplications.dm3.tab"
    conda: "envs/snake-cnv.yaml"
    shell: """wget "http://humanparalogy.gs.washington.edu/dm3/data/dm3genomicSuperDup.tab" -O {output}"""

#
# Repeat mask reference with RepeatMasker and TRF.
#

rule pad_repeats_by_36bp:
    input: bed="dm6.repeats_without_RNA.bed", chromosomes="dm6.chrom.sizes"
    output: "dm6.repeats_without_RNA.36bp_padded.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "bedtools slop -i {input.bed} -g {input.chromosomes} -b 36 | bedtools merge -i stdin -d 0 > {output}"

rule download_dm6_chromosome_sizes:
    output: "dm6.chrom.sizes"
    conda: "envs/snake-cnv.yaml"
    shell: """wget "http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/{output}" -O {output}"""

rule merge_repeatmasker_without_RNA_and_trf_outputs:
    input: "dm6.trf.bed.gz", "simpleRepeat.bed.gz", "windowmaskerSdust.bed.gz", "dm6.repeatmasker_without_RNA.bed.gz"
    output: "dm6.repeats_without_RNA.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | cut -f 1-3 | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d 0 > {output}"

rule merge_repeatmasker_and_trf_outputs:
    input: "dm6.trf.bed.gz", "dm6.repeatmasker.bed.gz"
    output: "dm6.repeats.bed"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | cut -f 1-3 | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d 0 > {output}"

rule get_repeatmasker_output_without_RNA:
    input: "dm6.repeatmasker.bed.gz"
    output: "dm6.repeatmasker_without_RNA.bed.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | grep -iv RNA | bgzip -c > {output}"

rule convert_repeatmasker_out_to_bed:
    input: "dm6.fa.out.gz"
    output: "dm6.repeatmasker.bed.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | $HOME/src/fasta_tools/out_to_bed.sh /dev/stdin | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d 0 -c 4 -o distinct | bgzip -c > {output}"

rule convert_windowmasker_to_bed:
    input: "windowmaskerSdust.txt.gz"
    output: "windowmaskerSdust.bed.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | cut -f 2-4 | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d 0 | bgzip -c > {output}"

rule convert_simpleRepeat_to_bed:
    input: "simpleRepeat.txt.gz"
    output: "simpleRepeat.bed.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} | cut -f 2-4 | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d 0 | bgzip -c > {output}"

rule download_windowmasker_output:
    output: "windowmaskerSdust.txt.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/{output} {output}"

rule download_simpleRepeat_output:
    output: "simpleRepeat.txt.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/{output} {output}"

rule download_trf_output:
    output: "dm6.trf.bed.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/{output} {output}"

rule download_repeatmasker_output:
    output: "dm6.fa.out.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/{output} {output}"

#
# Prepare reference for alignments and variant calling.
#

rule index_reference_for_gatk:
    input: "reference/dm6.fasta"
    output: "reference/dm6.dict"
    conda: "envs/snake-cnv.yaml"
    shell: "samtools dict {input} > {output}"

rule index_reference_for_bwa:
    input: "reference/dm6.fasta"
    output: "reference/dm6.fasta.bwt"
    conda: "envs/snake-cnv.yaml"
    shell: "bwa index {input}"

rule index_reference_for_samtools:
    input: "reference/dm6.fasta"
    output: "reference/dm6.fasta.fai"
    conda: "envs/snake-cnv.yaml"
    shell: "samtools faidx {input}"

rule uncompress_reference:
    input: "reference/dm6.fasta.gz"
    output: "reference/dm6.fasta"
    conda: "envs/snake-cnv.yaml"
    shell: "zcat {input} > {output}"

rule download_reference:
    output: "reference/dm6.fasta.gz"
    conda: "envs/snake-cnv.yaml"
    shell: "rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz {output}"
