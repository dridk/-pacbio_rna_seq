import pandas as pd 
import seaborn as sns
import sys
from glob import glob 
import re
import matplotlib.pyplot as plt

GENE_NAME = sys.argv[1]
OUTPUT_FILE = sys.argv[2]

alldf = []
for file in glob(f"{GENE_NAME}.*.hash.bed"):
    match = re.findall(r"(\w+)\.(\w+).+", file)[0]
    gene, barcode = match
    df = pd.read_csv(file, sep="\t")
    df["gene"] = gene
    df["barcode"] = barcode
    alldf.append(df)
    

df = pd.concat(alldf)
df

ackdf = df.query(f"gene == '{gene}'")[["hash", "gene","barcode"]]
keep_hash = ackdf["hash"].value_counts().reset_index()
keep_hash = keep_hash[keep_hash["hash"] > 50]

keep_hash = keep_hash["index"]
keep_hash

subdf = ackdf[ackdf["hash"].isin(keep_hash)]


subdf = subdf.groupby(["barcode", "hash"])["gene"].count().reset_index().pivot(index="barcode", columns="hash", values="gene").fillna(0)


sns.clustermap(subdf)
plt.title(gene)
plt.savefig(OUTPUT_FILE)