import pandas as pd
from datetime import datetime
import os

today = datetime.now()
outDir = today.strftime('%Y-%m-%d') + "_processed"

barcodes = pd.read_csv("barcodes.csv", sep=",", header=None, names=["barName", "barFasta"])

metadata_pd = pd.read_csv("metadata.txt", sep=",", header=None, names=["sampleName", "barName"])
metadata_names = metadata_pd['sampleName'].tolist()
metadata_bar = metadata_pd['barName'].tolist()

metadata_headers = metadata_pd['sampleName'] + "_" + metadata_pd['barName']

barcodeSubset = barcodes[barcodes['barName'].isin(metadata_bar)]

# Test if metadata barcodes are equal to what has been found in barcode file
if barcodeSubset["barName"].tolist()==metadata_bar:
    print("All barcodes found. Resuming")
elif barcodeSubset["barName"].tolist()!=metadata_bar:
    print("Some barcodes in metadata could not be found. Have you misspelled or copy-pasted wrong?")
    raise

barcodeUsedDF = pd.merge(barcodeSubset, metadata_pd, on='barName', how='inner')

barcodeUsedDF["sampleHeader"] = '>' + barcodeUsedDF['sampleName'].map(str) + '_' + barcodeUsedDF['barName'].map(str)


barcodeUsedOut = pd.DataFrame(barcodeUsedDF, columns = ["sampleHeader","barFasta"])

samplenames = barcodeUsedDF['sampleName'].map(str)+ '_' + barcodeUsedDF['barName'].map(str)

samplenames.to_csv(os.path.join(outDir, "samplenames.txt"), mode="w", header=False, index = False, sep="\n")

barcodeUsedOut.to_csv(os.path.join(outDir, "barcodes_used.fasta"), mode="w", header=False, index = False, sep="\n", line_terminator="\n")