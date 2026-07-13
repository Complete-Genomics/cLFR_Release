configfile: "consensus_fasta.yaml"

from pathlib import Path

SAMPLE_ID = config["sample_id"]
CHROMS = config["chroms"]
NUM_SPLITS = int(config["params"].get("num_splits", 5))
SPLITS = list(range(NUM_SPLITS))

wildcard_constraints:
    split_idx = r"\d+"

WORKFLOW_DIR = Path(workflow.basedir)
REPO_DIR = WORKFLOW_DIR.parent
SCRIPT_DIR = str(REPO_DIR / "scripts")
OUTPUT_DIR = config["paths"].get("output_dir", "Align")
KEEP_ALIGN_DIR = "keep/Align"
SPLIT_BAM_DIR = "Make_Vcf/step3_hapcut/step1_modify_bam"
TMP_DIR = f"{OUTPUT_DIR}/tmp"
CONSENSUS_DIR = f"{OUTPUT_DIR}/consensus"

SEQUENCE_TYPE = config["params"].get("sequence_type", "pe").lower()
MAPPER = config["mapping"].get("mapper", "star").lower()


rule all:
    input:
        f"{CONSENSUS_DIR}/consensus.fixRC.fasta",
        f"{CONSENSUS_DIR}/consensus_frag_length_distribution.pdf",


def input_read1(wildcards):
    if SEQUENCE_TYPE == "pe":
        return config["inputs"]["read1"]
    return []


rule stage_fastqs:
    input:
        read1=input_read1,
        read2=config["inputs"]["read2"],
    output:
        read1="data/read_1.fq.gz",
        read2="data/read_2.fq.gz",
    params:
        sequence_type=SEQUENCE_TYPE,
    shell:
        """
        mkdir -p data
        if [ "{params.sequence_type}" = "pe" ]; then
            ln -sf "$(realpath {input.read1})" {output.read1}
        else
            : > {output.read1}
        fi
        ln -sf "$(realpath {input.read2})" {output.read2}
        """


rule split_reads:
    input:
        read1="data/read_1.fq.gz",
        read2="data/read_2.fq.gz",
    output:
        read1="data/split_read.1.fq.gz",
        read2="data/split_read.2.fq.gz",
        log="split_stat_read1.log",
    threads:
        config["threads"].get("split", 8)
    params:
        python=config["tools"].get("python", "python"),
        read_len=config["params"]["read_len"],
        read_len_r1=config["params"].get("read_len_r1", 0),
        cbc_len=config["params"]["umi_len"],
        bc_len_redundant=config["params"].get("bc_len_redundant", config["params"]["umi_len"]),
        bc_start=config["params"].get("bc_start", 0),
        gdna_start=config["params"].get("gdna_start", 1),
        additional_bc_start=config["params"].get("additional_bc_start", 0),
        additional_bc_len=config["params"].get("additional_bc_len", 0),
        gdna_start_r1=config["params"].get("gdna_start_r1", 1),
        additional_bc_len_r1=config["params"].get("additional_bc_len_r1", 0),
        pigz=config["tools"].get("pigz", "pigz"),
        sequence_type=SEQUENCE_TYPE,
    shell:
        """
        mkdir -p data
        is_split=$({params.python} -c 'import gzip,re,sys; path=sys.argv[1]; umi_len=int(sys.argv[2]); opener=gzip.open if path.endswith((".gz", ".gzip")) else open; handle=opener(path, "rt", errors="replace"); header=handle.readline().strip(); handle.close(); print(1 if re.search(r"#[ACGTNacgtn]{{%d,}}(?:[#/\\s]|$)" % umi_len, header) else 0)' {input.read2} {params.cbc_len})
        if [ "$is_split" = "1" ]; then
            if [ "{params.sequence_type}" = "pe" ]; then
                ln -sf "$(realpath {input.read1})" {output.read1}
            else
                : > {output.read1}
            fi
            ln -sf "$(realpath {input.read2})" {output.read2}
            echo "Detected UMI in FASTQ read header; skipped split_barcode_LFR.py." > {output.log}
        else
            {params.python} {SCRIPT_DIR}/split_barcode_LFR.py \
                --r1 {input.read1} \
                --r2 {input.read2} \
                --read_len {params.read_len} \
                --output data/split_read \
                --cbc_len {params.cbc_len} \
                --bc_len_redundant {params.bc_len_redundant} \
                --bc_start {params.bc_start} \
                --gdna_start {params.gdna_start} \
                --additional_bc_start {params.additional_bc_start} \
                --additional_bc_len {params.additional_bc_len} \
                --gdna_start_r1 {params.gdna_start_r1} \
                --read_len_r1 {params.read_len_r1} \
                --additional_bc_len_r1 {params.additional_bc_len_r1} \
                --threads {threads} \
                --pigz {params.pigz} \
                2> data/split_stat_read.err
        fi
        """


rule trim_reads:
    input:
        read1="data/split_read.1.fq.gz",
        read2="data/split_read.2.fq.gz",
    output:
        read1="data/split_read_1_trimmed.fastq.gz",
        read2="data/split_read_2_trimmed.fastq.gz",
    params:
        bbduk=config["tools"].get("bbduk", "bbduk.sh"),
        sequence_type=SEQUENCE_TYPE,
    shell:
        """
        if [[ "{params.sequence_type}" == "pe" ]]; then
            {params.bbduk} in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} qtrim=rl
        elif [[ "{params.sequence_type}" == "se" ]]; then
            {params.bbduk} in={input.read2} out={output.read2} qtrim=rl && touch {output.read1}
        else
            echo "Unknown sequence_type: {params.sequence_type}" >&2
            exit 1
        fi
        """


def mapping_ref(wildcards):
    if MAPPER in ("bwa", "minimap2"):
        return config["mapping"].get("ref_fasta", config["paths"]["ref_fasta"])
    return []


def mapping_read1(wildcards):
    if MAPPER == "minimap2":
        return "data/split_read_1_trimmed.fastq.gz"
    return "data/split_read.1.fq.gz"


def mapping_read2(wildcards):
    if MAPPER == "minimap2":
        return "data/split_read_2_trimmed.fastq.gz"
    return "data/split_read.2.fq.gz"


rule map_reads:
    input:
        ref=mapping_ref,
        read1=mapping_read1,
        read2=mapping_read2,
    output:
        bam=f"{KEEP_ALIGN_DIR}/{SAMPLE_ID}.sort.bam",
    threads:
        config["threads"].get("mapping", 16)
    params:
        platform=config["params"].get("platform", "BGI-seq"),
        bwa=config["tools"].get("bwa", "bwa"),
        bwa_mem=config["mapping"].get("bwa_mem", "mem"),
        star=config["tools"].get("star", "STAR"),
        star_index=config["mapping"].get("star_index", ""),
        hisat2=config["tools"].get("hisat2", "hisat2"),
        hisat2_index=config["mapping"].get("hisat2_index", ""),
        hisat2_splicesites=config["mapping"].get("hisat2_splicesites", ""),
        minimap2=config["tools"].get("minimap2", "minimap2"),
        minimap2_preset=config["mapping"].get("minimap2_preset") or "splice:sr",
        minimap2_anno_bed=config["mapping"].get("minimap2_anno_bed", ""),
        minimap2_sort_mem=config["mapping"].get("minimap2_sort_mem", "2G"),
        samtools=config["tools"].get("samtools", "samtools"),
        mapper=MAPPER,
        sequence_type=SEQUENCE_TYPE,
        sample_id=SAMPLE_ID,
    shell:
        r"""
        mkdir -p {KEEP_ALIGN_DIR} {OUTPUT_DIR}
        if [ "{params.mapper}" = "star" ]; then
            if [ "{params.sequence_type}" = "pe" ]; then
                {params.star} --runThreadN {threads} \
                    --genomeDir {params.star_index} \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFileNamePrefix {OUTPUT_DIR}/ \
                    --outFilterScoreMinOverLread 0 \
                    --outFilterMatchNminOverLread 0 \
                    --outFilterMatchNmin 0 \
                    --readFilesCommand zcat \
                    --readFilesIn {input.read1} {input.read2}
            else
                {params.star} --runThreadN {threads} \
                    --genomeDir {params.star_index} \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFileNamePrefix {OUTPUT_DIR}/ \
                    --outFilterScoreMinOverLread 0 \
                    --outFilterMatchNminOverLread 0 \
                    --outFilterMatchNmin 0 \
                    --readFilesCommand zcat \
                    --readFilesIn {input.read2}
            fi
            mv {OUTPUT_DIR}/Aligned.sortedByCoord.out.bam {output.bam}
        elif [ "{params.mapper}" = "hisat2" ]; then
            if [ "{params.sequence_type}" = "pe" ]; then
                {params.hisat2} -x {params.hisat2_index} \
                    -1 {input.read1} -2 {input.read2} \
                    --known-splicesite-infile {params.hisat2_splicesites} \
                    --dta-cufflinks --score-min L,0,-0.2 \
                    --max-intronlen 1000000 --pen-noncansplice 30 \
                    --add-chrname --no-softclip \
                    -p {threads} \
                | {params.samtools} sort -o {output.bam} -@ {threads} -O bam -
            else
                {params.hisat2} -x {params.hisat2_index} \
                    -U {input.read2} \
                    --known-splicesite-infile {params.hisat2_splicesites} \
                    --dta-cufflinks --score-min L,0,-0.2 \
                    --max-intronlen 1000000 --pen-noncansplice 30 \
                    --add-chrname --no-softclip \
                    -p {threads} \
                | {params.samtools} sort -o {output.bam} -@ {threads} -O bam -
            fi
        elif [ "{params.mapper}" = "minimap2" ]; then
            SORT_TMP=/dev/shm/minimap_tmp_$$
            echo "minimap2 samtools sort tmp: $SORT_TMP" >&2
            echo "minimap2 preset: {params.minimap2_preset}" >&2
            mkdir -p $SORT_TMP
            bed_arg=""
            if [ -n "{params.minimap2_anno_bed}" ]; then
                bed_arg="--junc-bed {params.minimap2_anno_bed}"
            fi
            {params.minimap2} -ax {params.minimap2_preset} \
                -t {threads} --secondary=no --sam-hit-only \
                $bed_arg {input.ref} {input.read2} \
            | {params.samtools} view -@ 4 -b - \
            | {params.samtools} sort -@ {threads} -m {params.minimap2_sort_mem} \
                -T $SORT_TMP/sort \
                -o {output.bam} -
            rm -rf $SORT_TMP
        elif [ "{params.mapper}" = "bwa" ]; then
            rg="@RG\tID:{params.sample_id}\tSM:{params.sample_id}\tPL:{params.platform}"
            if [ "{params.sequence_type}" = "pe" ]; then
                {params.bwa} {params.bwa_mem} -M -R "$rg" -C \
                    -t {threads} {input.ref} {input.read1} {input.read2} \
                | {params.samtools} sort -o {output.bam} -@ {threads} -O bam -
            else
                {params.bwa} {params.bwa_mem} -M -R "$rg" -C \
                    -t {threads} {input.ref} {input.read2} \
                | {params.samtools} sort -o {output.bam} -@ {threads} -O bam -
            fi
        else
            echo "Unknown mapper: {params.mapper}" >&2
            exit 1
        fi
        """


rule mark_dups:
    input:
        bam=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.bam",
    output:
        bam=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.markdup.bam",
        metrics=f"{OUTPUT_DIR}/{{sample_id}}_dedup_metrics.txt",
    params:
        gatk=config["tools"].get("gatk", "gatk"),
        bc_condition=config["params"].get("bc_condition", "random_bc"),
        skip_markdup=config["params"].get("skip_markdup", False),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        if [[ "{params.skip_markdup}" == "True" ]]; then
            ln -sf "$(realpath {input.bam})" {output.bam}
            : > {output.metrics}
        elif [[ "{params.bc_condition}" == *"random_bc"* ]]; then
            {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -BARCODE_TAG BX -M {output.metrics}
        elif [[ "{params.bc_condition}" == *"standard"* ]]; then
            {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -BARCODE_TAG BC -M {output.metrics}
        else
            {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics}
        fi
        """


rule index_sort_bam:
    input:
        bam=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.bam",
    output:
        bai=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.bam.bai",
    threads:
        config["threads"].get("mapping", 16)
    params:
        samtools=config["tools"].get("samtools", "samtools"),
    shell:
        "{params.samtools} index -@ {threads} {input.bam} {output.bai}"


rule index_mark_dups:
    input:
        bam=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.markdup.bam",
    output:
        bai=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.markdup.bam.bai",
    threads:
        config["threads"].get("mapping", 16)
    params:
        samtools=config["tools"].get("samtools", "samtools"),
    shell:
        "{params.samtools} index -@ {threads} {input.bam} {output.bai}"


rule remove_duplicates:
    input:
        bam=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.markdup.bam",
        bai=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.markdup.bam.bai",
        prebam=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.bam",
        prebai=f"{KEEP_ALIGN_DIR}/{{sample_id}}.sort.bam.bai",
    output:
        bam=f"{OUTPUT_DIR}/{{sample_id}}.sort.removedup_rm000.bam",
        bai=f"{OUTPUT_DIR}/{{sample_id}}.sort.removedup_rm000.bam.bai",
    threads:
        config["threads"].get("mapping", 16)
    params:
        samtools=config["tools"].get("samtools", "samtools"),
        bc_condition=config["params"].get("bc_condition", "random_bc"),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        if [[ "{params.bc_condition}" == "random_bc" || "{params.bc_condition}" == "random_bc_umi_rc" || "{params.bc_condition}" == "pcrfree" ]]; then
            ln -sf "$(realpath {input.prebam})" {output.bam}
            ln -sf "$(realpath {input.prebai})" {output.bai}
        elif [[ "{params.bc_condition}" == *"standard"* ]]; then
            {params.samtools} view -@ {threads} -b -h -F 0x400 {input.bam} > {output.bam}
            {params.samtools} index -@ {threads} {output.bam} {output.bai}
        else
            ln -sf "$(realpath {input.bam})" {output.bam}
            ln -sf "$(realpath {input.bai})" {output.bai}
        fi
        """


rule get_chr_dict:
    input:
        bam=f"{OUTPUT_DIR}/{SAMPLE_ID}.sort.removedup_rm000.bam",
        bai=f"{OUTPUT_DIR}/{SAMPLE_ID}.sort.removedup_rm000.bam.bai",
    output:
        f"{OUTPUT_DIR}/samtools_idx.txt",
    params:
        samtools=config["tools"].get("samtools", "samtools"),
    shell:
        "{params.samtools} idxstats {input.bam} > {output}"


rule split_chrom_bam:
    input:
        bam=f"{OUTPUT_DIR}/{SAMPLE_ID}.sort.removedup_rm000.bam",
        bai=f"{OUTPUT_DIR}/{SAMPLE_ID}.sort.removedup_rm000.bam.bai",
    output:
        bam=f"{SPLIT_BAM_DIR}/{{sample_id}}_sort.markdup_{{chrom}}.bam",
        bai=f"{SPLIT_BAM_DIR}/{{sample_id}}_sort.markdup_{{chrom}}.bam.bai",
    params:
        samtools=config["tools"].get("samtools", "samtools"),
    shell:
        """
        mkdir -p {SPLIT_BAM_DIR}
        {params.samtools} view -bh {input.bam} {wildcards.chrom} > {output.bam}
        {params.samtools} index {output.bam} {output.bai}
        """


rule reformat_readid:
    input:
        bam=f"{SPLIT_BAM_DIR}/{{sample_id}}_sort.markdup_{{chrom}}.bam",
        bai=f"{SPLIT_BAM_DIR}/{{sample_id}}_sort.markdup_{{chrom}}.bam.bai",
    output:
        f"{TMP_DIR}/{{sample_id}}_{{chrom}}.name.bam",
    params:
        python=config["tools"].get("python", "python"),
        seq_type=SEQUENCE_TYPE,
        umi_len=config["params"]["umi_len"],
    shell:
        "{params.python} {SCRIPT_DIR}/consensus_fasta_supp.py "
        "--module reformat_readid "
        "--seq_type {params.seq_type} "
        "--umi_len {params.umi_len} "
        "--input_bam {input.bam} "
        "--output_bam {output}"


rule sort_reformatted_bam:
    input:
        f"{TMP_DIR}/{{sample_id}}_{{chrom}}.name.bam",
    output:
        f"{TMP_DIR}/{{sample_id}}_{{chrom}}.name.sort.bam",
    params:
        samtools=config["tools"].get("samtools", "samtools"),
    shell:
        "{params.samtools} sort -@ 2 -n -o {output} {input}"


rule get_consensus_fasta:
    input:
        bam=f"{TMP_DIR}/{{sample_id}}_{{chrom}}.name.sort.bam",
        chr_dict=f"{OUTPUT_DIR}/samtools_idx.txt",
    output:
        f"{TMP_DIR}/{{chrom}}/{{sample_id}}_{{chrom}}_{{split_idx}}.fasta",
    threads: 1
    params:
        python=config["tools"].get("python", "python"),
        ref=config["paths"]["ref_fasta"],
        min_reads=config["params"]["min_reads"],
        downsample_ratio=config["params"].get("downsample_ratio", 1.0),
        batch_id_arg=lambda wildcards: (
            f"--batch_id {config['params'].get('batch_id')}"
            if config["params"].get("batch_id")
            else ""
        ),
        samtools=config["tools"].get("samtools", "samtools"),
        stringtie=config["tools"].get("stringtie", "stringtie"),
        temp_dir=config["params"].get("temp_dir", "/dev/shm/consensus"),
        use_samtools_reference_arg=(
            "--use_samtools_reference"
            if config["params"].get("use_samtools_reference", False)
            else ""
        ),
    shell:
        "{params.python} {SCRIPT_DIR}/consensus_fasta.py "
        "--bam {input.bam} "
        "--ref_fasta {params.ref} "
        "--output_fasta {output} "
        "--chrom {wildcards.chrom} "
        "--dict_file {input.chr_dict} "
        "--split_index {wildcards.split_idx} "
        "--num_splits {NUM_SPLITS} "
        "--min_reads {params.min_reads} "
        "--downsample_ratio {params.downsample_ratio} "
        "{params.batch_id_arg} "
        "--samtools {params.samtools} "
        "--stringtie {params.stringtie} "
        "--temp_dir {params.temp_dir} "
        "{params.use_samtools_reference_arg}"


rule fix_consensus_format:
    input:
        f"{TMP_DIR}/{{chrom}}/{{sample_id}}_{{chrom}}_{{split_idx}}.fasta"
    output:
        f"{TMP_DIR}/{{chrom}}/{{sample_id}}_{{chrom}}_{{split_idx}}.noN.fix.fasta"
    params:
        python=config["tools"].get("python", "python")
    shell:
        "{params.python} {SCRIPT_DIR}/consensus_fasta_supp.py "
        "--module fix_fasta --input_fasta {input} --output_fasta {output}"


rule merge_consensus_fasta:
    input:
        expand(
            f"{TMP_DIR}/{{chrom}}/{{sample_id}}_{{chrom}}_{{split_idx}}.noN.fix.fasta",
            sample_id=SAMPLE_ID,
            chrom=CHROMS,
            split_idx=SPLITS,
        ),
    output:
        f"{CONSENSUS_DIR}/consensus.fasta",
    shell:
        "mkdir -p {CONSENSUS_DIR} && cat {input} > {output}"


rule map_fasta_consensus:
    input:
        f"{CONSENSUS_DIR}/consensus.fasta"
    output:
        f"{CONSENSUS_DIR}/consensus.paf"
    threads:
        config["threads"].get("mapping", 16)
    params:
        minimap2=config["tools"].get("minimap2", "minimap2"),
        ref=config["paths"].get("consensus_direction_ref_fasta", config["paths"]["ref_fasta"])
    shell:
        "{params.minimap2} -x asm20 -t {threads} {params.ref} {input} > {output}"


rule correct_direction_consensus:
    input:
        fasta=f"{CONSENSUS_DIR}/consensus.fasta",
        paf=f"{CONSENSUS_DIR}/consensus.paf"
    output:
        fasta=f"{CONSENSUS_DIR}/consensus.fixRC.fasta",
        not_mapped=f"{CONSENSUS_DIR}/frag_not_in_mapped.fixRC.fasta"
    params:
        python=config["tools"].get("python", "python"),
        correct_rc=config["params"].get("correct_rc", True),
        umi_len=config["params"].get("umi_len", 15),
        adapter_seq=config["params"].get("adapter_seq", ""),
        flank_end=config["params"].get("flank_end", 0)
    run:
        if params.correct_rc:
            shell(
                "{params.python} {SCRIPT_DIR}/fix_fasta_direction.py "
                "--input_fasta {input.fasta} "
                "--paf {input.paf} "
                "--output_fasta {output.fasta} "
                "--not_mapped_fasta {output.not_mapped} "
                "--umi_len {params.umi_len} "
                "--adapter_seq '{params.adapter_seq}' "
                "--flank_end {params.flank_end}"
            )
        else:
            shell("cp {input.fasta} {output.fasta} && : > {output.not_mapped}")


rule fasta_frag_len_distribution_consensus:
    input:
        f"{CONSENSUS_DIR}/consensus.fixRC.fasta",
    output:
        f"{CONSENSUS_DIR}/consensus_frag_length_distribution.pdf",
    params:
        python=config["tools"].get("python", "python"),
        minreads_fasta=config["params"]["minreads_fasta"],
    shell:
        "{params.python} {SCRIPT_DIR}/exon2fasta.py "
        "--fasta {input} "
        "--outdir {CONSENSUS_DIR} "
        "--name consensus "
        "--minreads_fasta {params.minreads_fasta} "
        "--module metrics_basic"
