#!/bin/bash
fastqDir=$1
threadNo="${2:-30}"
OUT_DIR=$(date '+%Y-%m-%d'_processed)
version="2.0"

set -o errexit -o pipefail -o noclobber -o nounset
[  -z "$fastqDir" ] && { echo "No path to a fastq directory provided. Exiting..."; exit 0; }

mkdir -p $OUT_DIR

python /opt/scripts/barcode_compiler.py

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

echo "Running your snakemake workflow ($workFlowSelect) with $threadNo threads"
snakemake --cores $threadNo --config fastq=$fastqDir threads=$threadNo workflow=$workFlowSelect
