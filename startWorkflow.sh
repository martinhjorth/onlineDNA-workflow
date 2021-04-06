#!/bin/bash
#conda install -c conda-forge biopython
fastqDir=$1
threadNo="${2:-30}"
OUT_DIR=$(date '+%Y-%m-%d'_processed)
version="2.0"

set -o errexit -o pipefail -o noclobber -o nounset
[  -z "$fastqDir" ] && { echo "No path to a fastq directory provided. Exiting..."; exit 0; }


mkdir -p $OUT_DIR
# Clean module space
module purge

# Load Python
module load Python/3.6.4-foss-2018a --quiet

# Run barcode_compiler
python scripts/barcode_compiler.py

# Load snakemake environment
echo "Loading Miniconda"
module load Miniconda3/4.9.2-foss-2020b --quiet

# echo "Activating snakemake"
#conda init bash
set +eu
. activate /srv/MA/users/mha/CondaEnvs/snakemake
set -eu
#conda activate /srv/MA/users/mha/CondaEnvs/snakemake

# Select workflow
title="onlineDNA showcase workflow v. $version"
prompt="Choose which primer-set you have used:"
options=("bV1-3", "bV4-A")
echo "$title"
PS3="$prompt "
select opt in "${options[@]}" "Quit"; do 

    case "$REPLY" in
        1)
            echo ""
            echo "You chose bV1-3. Trimming reads accordingly"
            workFlowSelect=onlineDNA_V13_v2.0;
            break
            ;;
        2)
            echo ""
            echo "You chose bV4-A. Trimming reads accordingly"
            workFlowSelect=onlineDNA_V4_v2.0;
            break
            ;;
            $(( ${#options[@]}+1 )) ) echo "Goodbye!"; break;;
    *) echo "Invalid option. Try another one.";continue;;
  esac
done



#conda activate snakemake
echo "Running your snakemake workflow ($workFlowSelect) with $threadNo threads"
snakemake --cores 20 --config fastq=$fastqDir threads=$threadNo workflow=$workFlowSelect --use-conda