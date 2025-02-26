rule all:
    input:
        "results/non-target_homologous.csv"  # Это правило указывает, что конечным результатом всего пайплайна является файл 'non-target_homologous.csv'
        
rule download_genome:
    output:
        "data/GCF_000001405.13_GRCh37_genomic.fna.gz"  # Скачивание генома в формате .gz
    conda:
        "envs/genome_env.yaml"  # Используется conda окружение для этой задачи
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz -O {output}
        """  # Скачивание генома с ftp-ресурса

rule unzip_genome:
    input:
        "data/GCF_000001405.13_GRCh37_genomic.fna.gz"  # Исходный файл в формате .gz
    output:
        "data/GCF_000001405.13_GRCh37_genomic.fna"  # Разархивированный геном
    conda:
        "envs/unzip_env.yaml"  # Используется conda окружение для разархивирования
    shell:
        """
        gunzip -c {input} > {output}
        """  # Разархивирование файла .gz с помощью gunzip

rule extract_sequences:
    input:
        genome="data/GCF_000001405.13_GRCh37_genomic.fna",  # Геном в формате .fna
        bed="data/IAD143293_241_Designed_formatRefSeq.bed"  # Файл с диапазонами для извлечения последовательностей
    output:
        "results/target_regions_RefSeq.fasta"  # Извлеченные последовательности, сохраненные в формате .fasta
    conda:
        "envs/bedtools_env.yaml"  # Используется conda окружение для работы с bedtools
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output}
        """  # Извлечение последовательностей из генома по указанным в файле .bed регионам

rule create_blast_db:
    input:
        "data/GCF_000001405.13_GRCh37_genomic.fna"  # Исходный геномный файл для создания базы данных BLAST
    output:
        db_files=expand("results/my_blast_db.{ext}", ext=["nhr", "nin", "nsq", "ndb", "nto", "ntf", "not", "njs"])  # Набор файлов базы данных BLAST
    conda:
        "envs/blast_env.yaml"  # Используется conda окружение для создания базы данных BLAST
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out results/my_blast_db
        """  # Создание базы данных BLAST для поиска по нуклеотидным последовательностям

rule run_blast:
    input:
        db_files=expand("results/my_blast_db.{ext}", ext=["nhr", "nin", "nsq", "ndb", "nto", "ntf", "not", "njs"]),  # Все файлы базы данных BLAST
        query="results/target_regions_RefSeq.fasta"  # Запрос (извлеченные последовательности)
    output:
        "results/blast_results_RefSeq.txt"  # Результаты поиска BLAST
    conda:
        "envs/blastn_env.yaml"  # Используется conda окружение для выполнения поиска BLAST
    shell:
        """ 
        blastn -db results/my_blast_db -query {input.query} -out {output} -outfmt 6
        """  # Выполнение поиска BLAST с параметром outfmt 6 для текстового вывода результатов

rule filter_blast_results:
    input:
        "results/blast_results_RefSeq.txt"  # Результаты BLAST
    output:
        "results/regions_100Ident.csv"  # Фильтрация результатов для сохранения только тех, что с идентичностью 100%
    shell:
        """
        awk -F'\t' '$3 == "100.000" && $4 == ($10 - $9 + 1) && $7 == 1' results/blast_results_RefSeq.txt > results/regions_100Ident.csv
        """  # Фильтрация результатов BLAST для выбора регионов с полным совпадением

rule filter_non_target_homologous:
    input:
        "results/regions_100Ident.csv"  # Файл с регионами с 100% идентичностью
    output:
        "results/non-target_homologous.csv"  # Фильтрация, чтобы найти те регионы, которые являются гомологами, но не целевыми
    shell:
        """
        awk '{{print $1}}' {input} | sort | uniq -d | while read val; do awk -v var="$val" '$1 == var' {input}; done > {output}
        """  # Извлечение дублирующихся значений из первого столбца, затем поиск соответствующих строк в исходном файле
