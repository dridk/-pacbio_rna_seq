import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

import re 

# ACKR1.transcript.bc1001.unique.bed
INPUT_FILE = sys.argv[1]
LIMIT  = int(sys.argv[2])




name, barcode = re.findall(r"(\w+)\.transcript\.(bc\d+)\.unique\.bed", INPUT_FILE)[0]

OUTPUT_FILE = f"{name}.transcript.{barcode}.dist.png"


try:
	df = pd.read_csv(INPUT_FILE, sep="\t", header=None).iloc[:,0:4]
	df.columns = ["chr","start","end","count"]
	df = df.reset_index()
	df["freq"] = (df["count"] / df["count"].sum() * 100).round(2)

	sns.barplot(x="index",y="freq",data=df.head(LIMIT))
	plt.title(f"{name} - {barcode}")
	plt.xlabel("unique transcripts")
	plt.ylabel("Frequence (%)")
	
	plt.savefig(OUTPUT_FILE)

except:
	plt.Figure()
	plt.savefig(OUTPUT_FILE)