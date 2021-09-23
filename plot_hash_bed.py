import matplotlib.pyplot as plt
import re 
import seaborn as sns 
import pandas as pd 
import sys

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]

match = re.search(r"(\w+).(\w+).hash.bed", INPUT_FILE)

title = match.group(1) + " " + match.group(2)

try:
	subdf = pd.read_csv(INPUT_FILE, sep="\t")
	pldf = subdf["hash"].value_counts().reset_index()
	pldf["color"] = "#"+pldf["index"]
	pldf = pldf[pldf["hash"] > 1]
	sns.barplot(x="index",y="hash",data=pldf, palette=pldf["color"])
	plt.title(title)
	plt.xticks(rotation=90)
	plt.tight_layout()
	plt.savefig(OUTPUT_FILE)

except:
	plt.Figure()
	plt.savefig(OUTPUT_FILE)
