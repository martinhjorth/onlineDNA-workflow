import argparse
from Bio import SeqIO, SeqRecord
#from biopython import SeqIO, SeqRecord
import os
from typing import List
parser = argparse.ArgumentParser(description='Reverse complementing!')
parser.add_argument("-o", "--outputlocation", help="output files to this location", action="store")
parser.add_argument("-i", "--inputlocation", help="output files to this location", action="store")
parser.add_argument("-m", "--metadatafile", help="output file name and path eg. ./output.csv", action="store")
args = parser.parse_args()
output_directory = args.outputlocation
input_directory = args.inputlocation
metadata_file_location = args.metadatafile

#input_directory = snakemake.input[0]
#metadata_file_location = snakemake.input[1]
#output_directory = snakemake.output[0]

class Metadata:
    library_id: str
    def __init__(self, name: str, barcode: str):
        self.name = name
        self.barcode = barcode
    def set_library_id(self, name: str):
        self.library_id = name
def read_metadata_file(path: str) -> List[Metadata]:
    text = open(path, "r").read().splitlines()
    metadata_items: List[Metadata] = []
    for line in text:
        parts = line.split(',')
        metadata_items.append(Metadata(parts[0], parts[1]))
    return metadata_items
def find_temp_file_with_extension_and_suffix(path: str, extension: str, suffix: str) -> tuple:
    for root, directories, files in os.walk(path, followlinks=True):
        for file in files:
            if file.endswith(suffix + "." + extension):
                return root, file
def return_string_proceeding_match(string: str, match: str) -> str:
    return string.rpartition(match)[0]
def record_should_be_reversed(record: SeqRecord) -> bool:
    return record.description.endswith(" rc")
def get_time_from_description(record: SeqRecord) -> str:
    return [i for i in record.description.split(" ") if i.startswith('start_time=')][0].split("=")[1]
def output_fastq_records_to_file(records: [SeqRecord], path: str):
    with open(path, "w") as output_handle:
        SeqIO.write(records, output_handle, "fastq")
def create_custom_description(record: SeqRecord) -> str:
    return "time=" + get_time_from_description(record) + "barcode=" + metadata_pair.barcode
def create_output_directory():
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
metadata_file_records = read_metadata_file(metadata_file_location)
create_output_directory()
for metadata_pair in metadata_file_records:
    metadata_temp_file_path = find_temp_file_with_extension_and_suffix(input_directory, "tmp", metadata_pair.barcode)
    if metadata_temp_file_path:
        metadata_pair.set_library_id(return_string_proceeding_match(metadata_temp_file_path[1], "_UDP"))
        processed_records = []
        tmp_file_path = os.path.join(metadata_temp_file_path[0], metadata_temp_file_path[1])
        for record in SeqIO.parse(tmp_file_path, "fastq"):
            if record_should_be_reversed(record):
                record = record.reverse_complement(True, True, True, True, True, True, True)
            record.description = create_custom_description(record)
            record.id = (record.id + "_" + record.description) # SeqIO by default makes everything up until first space the "id" and the rest of the header is description. Here we have trimmed the header and put it back together as "id"
            processed_records.append(record)
        output_fastq_records_to_file(processed_records, os.path.join(output_directory, metadata_pair.name + ".fastq"))
