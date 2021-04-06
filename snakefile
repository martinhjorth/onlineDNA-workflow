import os, glob, subprocess, shutil, Bio#, argparse, random, re, gzip, shutil, subprocess
from datetime import datetime
#from pick import pick
import pandas as pd

dbFasta = "/srv/MA/users/mha/databases/midas37/midas37_notax.fa" # Path to database fasta
dbTax = "/srv/MA/users/mha/databases/midas37/tax_complete_qiime.txt" # # Path to database taxonomy
dbName = dbFasta.split("/")[-1].split("_")[-2]

today = datetime.now().strftime('%Y-%m-%d')
outDir = today + "_processed"

# Get input variables
fastqDir = config['fastq']
threadNo = int(config['threads'])
workflowSelected = config['workflow']

with open(os.path.join(outDir, "samplenames.txt"), "r") as f:
    sampleNames = f.readlines()
samples_list = [x.strip() for x in sampleNames]
libnames = [x.split("_")[0] for x in samples_list]

rule all:
    input:
        expand([os.path.join(outDir, "reads_trim.fq"),
                os.path.join(outDir, "cutadapt_outAdapt_qc.txt"),
                os.path.join(outDir, "demultiplexed/{libID}.fastq"),
                os.path.join(outDir, "mapped/{libID}.sam"),
                os.path.join(outDir, "mapped/{libID}.idmapped.txt"),
                os.path.join(outDir, today + "_mapping-table_" + dbName + ".txt"),
                os.path.join(outDir, "mapped/temp/all_mappings_filt25.txt"),
                os.path.join(outDir, "mapped/temp/new_mappings.txt"),
                os.path.join(outDir, "mapped/temp/new_mappings_filt25.txt"),
                os.path.join(outDir, today + "_stats.txt")], name=samples_list, libID=libnames)

rule cat_fastq:
    input: fastqDir
    output: temp(os.path.join(outDir, "tempReads.fq"))
    run:
        with open(output[0], 'wb') as outfile:
            for filename in glob.iglob(os.path.join(input[0], '*.fastq')):
                try:
                    with open(filename, 'rb') as infile:
                        shutil.copyfileobj(infile, outfile)
                except IOError:
                    raise

# Set trimming parameters based on selected workflow
if workflowSelected == "onlineDNA_V4_v2.0":
    ampMin_len = 275
    ampMax_len = 450
elif workflowSelected == "onlineDNA_V13_v2.0":
    ampMin_len = 500
    ampMax_len = 800


rule trim_outerAdapt:
    """ Trim end adaptors and filter by length.
     -- 20 bp of each terminal adaptor is targeted
     -- Outer adaptors are searched for together with cutadapt <ADP1_outer>...<ADP2_outer>
     -- Sequences are reverse complemented to obtain same i7i5 adapter orientation """
    input:
        os.path.join(outDir, "tempReads.fq")
    output:
        fastq=os.path.join(outDir, "reads_trim.fq"),
        qc=os.path.join(outDir, "cutadapt_outAdapt_qc.txt")
    params:
        ampMax_len = ampMax_len,
        ampMin_len = ampMin_len,
        threads = threadNo
    message: "Trimming outer adapters (Illumina pads)"
    log: os.path.join(outDir, "logs/cutadapt/outer_adapters.log")
    shell:
        """
        #!/bin/bash
        module load cutadapt/2.8-foss-2018a-Python-3.6.4 --quiet
        cutadapt -j {params.threads} -e 0.20 -O 10 -m {params.ampMin_len} -M {params.ampMax_len} --discard-untrimmed --revcomp -g CAGAAGACGGCATACGAGAT...GTGTAGATCTCGGTGGTCGC -o {output.fastq} {input} > {output.qc}
        """

rule demultiplex:
    """ Demultiplex Nextera-styled barcodes using the compiled barcode fasta file
    -- Dual barcode demultiplexing requires i7i5 adaptor orientation
    -- If orientation is not streamlined reads might be wrongly assigned as
    -- some dual barcoding schemes re-use barcode sequences in reverse complement orientation """
    input:
        barcodeFile=os.path.join(outDir, "barcodes_used.fasta"),
        fastqFile=os.path.join(outDir, "reads_trim.fq")
    output: temp(expand(os.path.join(outDir, "{name}.tmp"), name=samples_list))
    params:
        ampMax_len = ampMax_len,
        ampMin_len = ampMin_len,
        outDir = outDir,
        threads = threadNo
    message: "Demultiplexing trimmed reads"
    log: "logs/cutadapt/demultiplexing.log"
    shell:
        """
        #!/bin/bash
        module load cutadapt/2.8-foss-2018a-Python-3.6.4 --quiet
        cutadapt -e 0.20 -O 10 -m {params.ampMin_len} -M {params.ampMax_len} -g file:{input.barcodeFile} -o {params.outDir}/{{name}}.tmp {input.fastqFile}
        """

rule revcomp_trim:
    """ Trim internal adapter sequences and revert reverse complement
    -- Trim of inner adaptors based on adaptor length
    -- Revert sequence orientation to prepare for downstream consensus calling """
    input:
        demux_files = expand(os.path.join(outDir, "{name}.tmp"), name=samples_list),
        script = "scripts/revcomp.py",
        demuxDir = outDir,
        metadata="metadata.txt"
    output: os.path.join(outDir, "demultiplexed/{libID}.fastq")
    params:
        demuxDir = lambda wildcards: os.path.join(outDir, "demultiplexed")
    message: "Trimming internal adapters and reverse complementing"
    shell:
        """
        #module load Python/3.6.4-foss-2018a
        python {input.script} -i {input.demuxDir} -m {input.metadata} -o {params.demuxDir}
        """

rule map_files:
    """ Map processed fastq files and filter SAM flags 4, 256 and 2048 away (unmapped, non-primary and supplementary)"""
    input:
        fastqFiles = os.path.join(outDir, "demultiplexed/{libID}.fastq"),
        dbFasta_input = dbFasta
    output: os.path.join(outDir, "mapped/{libID}.sam")
    params: threads = threadNo
    message: "Mapping sample {wildcards.libID}"
    shell:
        """
        module load Minimap2/2.17-foss-2018a --quiet
        module load SAMtools/1.10-foss-2018a --quiet
        minimap2 -ax map-ont -t {params.threads} --secondary=no {input.dbFasta_input} {input.fastqFiles} |
        samtools view -F 4 -F 256 -F 2048 - -o {output}
        """

rule mapping_essentials:
    """ Extract fields of interest from SAM files and output in txt """
    input: os.path.join(outDir, "mapped/{libID}.sam")
    output: os.path.join(outDir, "mapped/{libID}.idmapped.txt")
    message: "Extracting essential information from mapped SAM files"
    shell:
        """
        #!/bin/bash
        sed '/^@/ d' {input} |
        awk '{{
        for(i=1;i<=NF;i++){{
        if($i ~ /^NM:i:/){{sub("NM:i:", "", $i); mm = $i}}
        }}
        split($6, count, /[^0-9]+/);
        split($6, type, /[^A-Z]*/);
        for(i=1; i<= length(count)-1; i++){{
        if(type[i + 1] ~ /[DIM]/){{aln+=count[i]}};
        }}
        print $1, $2, $3, length($10), aln, (aln - mm)/aln, $12, $14, $20
        aln=0;
        }}' |
        sed 's/_time=/\\t/' | sed 's/Zbarcode=/\\t/' > {output}
        """
rule make_map_table:
    """ Load the mappings into R, merge with taxonomy and output a mapping table """
    input:
        dbTax = dbTax,
        newMappings_dir = outDir,
        newMappings = expand(os.path.join(outDir, "mapped/{libID}.idmapped.txt"), libID=libnames)
    output:
        maptable = os.path.join(outDir, today + "_mapping-table_" + dbName + ".txt"),
        all_mappings = os.path.join(outDir, "mapped/temp/all_mappings_filt25.txt"),
        new_mappings = os.path.join(outDir, "mapped/temp/new_mappings.txt"),
        new_mappingsFilt = os.path.join(outDir, "mapped/temp/new_mappings_filt25.txt")
    params:
        dbname = dbName
    message:
        "Creating mapping table of current and previous sequencing runs"
    script:
        "scripts/create_mapping_tables.R"

rule qc_report:
    input:
        newMappings = os.path.join(outDir, "mapped/temp/new_mappings_filt25.txt")
    output: os.path.join(outDir, today + "_stats.txt")
    run:
        import pandas as pd
        mapTable = pd.read_csv(input.newMappings, sep="\t", header=0)
        readCounts = mapTable["SeqID"].value_counts().sort_index()
        refCounts = mapTable.groupby('SeqID')["OTU"].nunique()
        seqStats = pd.concat([readCounts, refCounts], axis=1).reset_index()
        seqStats.columns = ["SeqID", "ReadCount", "UniqueASVs"]
        seqStats.to_csv(output[0], mode="w", header=True, index = False, sep="\t")
