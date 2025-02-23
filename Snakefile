rule all:
    input:
        "results/homologous_regions_2.csv"

rule download_genome:
    output:
        "data/GCF_000001405.13_GRCh37_genomic.fna.gz"
    conda:
        "envs/genome_env.yaml"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz -O {output}
        """

rule unzip_genome:
    input:
        "data/GCF_000001405.13_GRCh37_genomic.fna.gz"
    output:
        "data/GCF_000001405.13_GRCh37_genomic.fna"
    conda:
        "envs/unzip_env.yaml"
    shell:
        """
        gunzip -c {input} > {output}
        """

rule extract_sequences:
    input:
        genome="data/GCF_000001405.13_GRCh37_genomic.fna",
        bed="data/IAD143293_241_Designed_formatRefSeq.bed"
    output:
        "results/target_regions_RefSeq.fasta"
    conda:
        "envs/bedtools_env.yaml"
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output}
        """

rule create_blast_db:
    input:
        "data/GCF_000001405.13_GRCh37_genomic.fna"
    output:
        db_files=expand("results/my_blast_db.{ext}", ext=["nhr", "nin", "nsq", "ndb", "nto", "ntf", "not", "njs"])
    conda:
        "envs/blast_env.yaml"
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out results/my_blast_db
        """

rule run_blast:
    input:
        db_files=expand("results/my_blast_db.{ext}", ext=["nhr", "nin", "nsq", "ndb", "nto", "ntf", "not", "njs"]),
        query="results/target_regions_RefSeq.fasta"
    output:
        "results/blast_results_RefSeq.txt"
    conda:
        "envs/blastn_env.yaml"
    shell:
        """ 
        blastn -db results/my_blast_db -query {input.query} -out {output} -outfmt 6
        """

rule filter_blast_results:
    input:
        "results/blast_results_RefSeq.txt"
    output:
        "results/homologous_regions_2.csv"
    shell:
        """
        awk -F'\t' '$3 == "100.000" && $4 == ($10 - $9 + 1)' results/blast_results_RefSeq.txt > results/homologous_regions_2.csv
        """

