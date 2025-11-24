import glob 
from os.path import basename,dirname
import os
from datetime import datetime

configfile: "aufame_config.yaml"

GENOMES_DIR = config["genomes_dir"]
ANNOTOOLS = config["annotools"]

EMAPPER2GBK_ENV = config["emapper2gbk_env"]
BAKTA_ENV = config["bakta_env"]
PROKKA_ENV = config["prokka_env"]
EGGNOG_ENV = config["eggnog_env"]
METAGE2METABO = config["metage2metabo_env"]

TAXFILE = config["taxfile"] 

EGGNOG_DB_PATH = "/dev/shm/eggnog_db_" + datetime.today().strftime("%Y%m%d_%H%M%S") if config["eggnog_db_optim_path"] == "" else config["eggnog_db_optim_path"]

workdir: config["output_dir"]

SAMPLES = [basename(dirname(file)) for file in glob.glob(f"{GENOMES_DIR}/*/*.fasta")]
INPUT_FASTA = lambda wildcards: f"{GENOMES_DIR}/{wildcards.sample}/{wildcards.sample}.fasta"

def create_taxon_file(wildcards, samples, taxfile):
    """
        From taxon file, generate another version of taxon file 
        interpretable for mpwt
        Input : 
            wildcards (str) : string corresponding to annotation tool (smk wildcards)
            samples (list) : samples names list
        Output : 
            taxfile written in annotool's directory
    """
    import pandas as pd

    annotool = wildcards
    df_taxons = pd.read_csv(taxfile, sep='\t')
    df_to_write = pd.DataFrame(columns = ["species", "taxon_id"])

    ## fill dfs
    col_filename = df_taxons.columns[2]
    col_taxID = df_taxons.columns[1]
    for index, row in df_taxons.iterrows() : 
        filename = row[col_filename]
        taxID = row[col_taxID]
        if filename in samples :
            df_to_write.loc[len(df_to_write)] = [filename,taxID]

    ## writing new file in annotation subdir
    tax_file = os.path.join(annotool, 'taxon_id.tsv')
    df_to_write.to_csv(tax_file, sep="\t", index=False)

rule all:
    input:
        "tsv_files/reactions.tsv"

rule prokka:
    input: 
        INPUT_FASTA
    output: 
        gbk = "prokka/{sample}/{sample}.gbk"
    conda: 
        PROKKA_ENV
    params:
        EXTENSIONS_TO_REMOVE = " ".join(config["to_remove"]["prokka"])
    shell: """
           mkdir -p prokka/{wildcards.sample}

           LANG=en_US.UTF-8
           prokka {input} --outdir prokka/{wildcards.sample} \
           --prefix {wildcards.sample} --compliant --force \
           > prokka/{wildcards.sample}/{wildcards.sample}_log.log
           
           for ext in {params.EXTENSIONS_TO_REMOVE}; do
                rm -f prokka/{wildcards.sample}/{wildcards.sample}$ext || true
           done
           """

rule bakta:
    input: 
        INPUT_FASTA
    output: 
        gbk = "bakta/{sample}/{sample}.gbk",
        faa = "bakta/{sample}/{sample}.faa"
    conda: 
        BAKTA_ENV
    params:
        BAK_DB=config["bak_db"],
        EXTENSIONS_TO_REMOVE = " ".join(config["to_remove"]["bakta"])
    resources:
        process_data_jobs = 1
    shell: """
           mkdir -p bakta/{wildcards.sample}

           bakta {input} --output bakta/{wildcards.sample} \
           --prefix {wildcards.sample} --db {params.BAK_DB} \
           --compliant --force
           
           mv bakta/{wildcards.sample}/{wildcards.sample}.gbff bakta/{wildcards.sample}/{wildcards.sample}.gbk
           
           for ext in {params.EXTENSIONS_TO_REMOVE}; do
                rm -f bakta/{wildcards.sample}/{wildcards.sample}$ext || true
           done
           """
        
rule eggnog:
    input:
        faa_bakta = "bakta/{sample}/{sample}.faa"
    output: 
        annot = temp("eggnog/{sample}/{sample}.emapper.annotations"),
        hits = temp("eggnog/{sample}/{sample}.emapper.hits"),
        ortho = temp("eggnog/{sample}/{sample}.emapper.seed_orthologs")
    conda: 
        EGGNOG_ENV
    threads: 
        8
    params:
        FILEBASE="eggnog/{sample}",
        EGG_DB=config["egg_db"],
        db_in_mem_path=EGGNOG_DB_PATH,
        db_in_mem_bool=config["eggnog_db_optim"],
        mem_threshold=50
    shell: """
        mem_threshold={params.mem_threshold}
        db_in_mem={params.db_in_mem_path}
        mkdir -p $db_in_mem
        tmp_dir="/projet/tmp/aufame"

        mkdir -p eggnog/{wildcards.sample}
        if [ {params.db_in_mem_bool} == "yes" ]; then
            if [ -z "$( ls -A $db_in_mem )" ]; then
                echo "WARNING : adding database to memory can raise permissions errors."
                echo "If so, please transfer yourself the database with the following command : "
                echo "cp {params.EGG_DB}/* $db_in_mem/"
                memfree=$(awk '/MemFree/ {{ printf "%.3f \\n", $2/1024/1024 }}' /proc/meminfo)
                memfree_ok=$(echo $memfree | awk -v thresh="$mem_threshold" '{{ print ($1 > thresh) ? "true" : "false" }}')
                if [ $memfree_ok == "true" ]; then
                    echo "Enough space for copying Eggnog DB in memory, processing..."
                    cp {params.EGG_DB}/* $db_in_mem/
                    options_db=$db_in_mem
                    echo "Launching Eggnog-mapper accordingly."
                else 
                    echo "Not enough space in memory, processing with '--dbmem' option."
                    options_db="{params.EGG_DB} --dbmem"
                fi
            else 
                echo "Eggnog DB in memory, processing accordingly..."
                options_db=$db_in_mem
            fi
        else 
            echo "No request of shifting Eggnog-DB to memory, processing it with '--dbmem' option."
            options_db="{params.EGG_DB} --dbmem"
        fi

        mkdir -p $tmp_dir

        emapper.py -i {input.faa_bakta} -o {wildcards.sample} \
            --itype proteins \
            --data_dir $options_db \
            --output_dir {params.FILEBASE} \
            --override \
            --cpu {threads}  --scratch_dir $tmp_dir --target_taxa 2 \
            > eggnog/{wildcards.sample}/{wildcards.sample}_log.log
        """

rule emapper2gbk:
    input:
        fasta = INPUT_FASTA, 
        fastap = "bakta/{sample}/{sample}.faa",
        annot = "eggnog/{sample}/{sample}.emapper.annotations"
    output: 
        gbk = "eggnog/{sample}/{sample}.gbk"
    conda: 
        EMAPPER2GBK_ENV
    shell: """
           emapper2gbk genes \
           -fn {input.fasta} \
           -fp {input.fastap} \
           -a {input.annot} \
           -o {output.gbk}
           """

rule create_taxon:
    input:
        TAXFILE
    output:
        # os.path.abspath(f"{wildcards.annotool}/taxon_id.tsv")
        os.path.join("{annotool}", "taxon_id.tsv")
    run:
        create_taxon_file(wildcards.annotool, SAMPLES, TAXFILE)

## expand concerns all of the samples to be created
## annotool is the only wildcard
rule mpwt:
    input:
        gbks = expand("{{annotool}}/{sample}/{sample}.gbk", sample=SAMPLES),
        taxon_file = os.path.join("{annotool}", "taxon_id.tsv")
    output:
        expand("mpwt/{{annotool}}/{sample}.zip", sample=SAMPLES)
    conda: 
        METAGE2METABO
    threads: 
        32
    resources:
        process_data_jobs = 6
    shell: """
        if [ {wildcards.annotool} == "eggnog" ]; then 
            rm -rf {EGGNOG_DB_PATH} || true
        fi 

        mkdir -p mpwt

        mpwt -f {wildcards.annotool}/ \
        --cpu {threads} \
        -o mpwt/{wildcards.annotool} \
        --patho --flat --clean --md -r -v
        """

rule pgdb2padmet:
    input: 
        "mpwt/{annotool}/{sample}.zip"
    output: 
        "padmet/{sample}/{sample}_{annotool}.padmet"
    conda: 
        METAGE2METABO
    params: 
        METACYC_REF=config["metacyc_ref"]
    shell: """
           mkdir -p padmet/

           unzip {input} -d mpwt/{wildcards.annotool}/{wildcards.sample}

           sleep 1m
           padmet pgdb_to_padmet --source=annot_{wildcards.annotool} \
           --pgdb=mpwt/{wildcards.annotool}/{wildcards.sample} \
           --output={output} --extract-gene --no-orphan \
           --padmetRef={params.METACYC_REF} -v
           """ 


rule merge_padmet:
    ## here, as annotools is not a wildcard in output, 
    ## we have to use an expand to assess all possibilities
    input:   
        expand("padmet/{{sample}}/{{sample}}_{annotool}.padmet", annotool=ANNOTOOLS)
    output: 
        "merged_padmet/{sample}.padmet"
    conda: 
        METAGE2METABO
    shell: """
           mkdir -p merged_padmet/
           
           sleep 1m
           padmet padmet_to_padmet --to_add=padmet/{wildcards.sample} --output={output} -v
           """


rule compare_padmet:
    input: 
        expand("merged_padmet/{sample}.padmet", sample=SAMPLES)
    output: 
        "tsv_files/reactions.tsv"
    conda: 
        METAGE2METABO
    shell: """
           mkdir -p tsv_files

           sleep 1m
           padmet compare_padmet --padmet=merged_padmet --output=tsv_files -v
           touch {output}
           """ 