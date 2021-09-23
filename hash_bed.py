import hashlib
import pandas as pd 
import sys

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]

df = pd.read_csv(INPUT_FILE, sep="\t", names = ["chr","start","end","name","qual","sens","a","b","col","count","estart","eend"])
subdf = df[["chr","start","end","count","estart","eend"]].copy()
subdf["hash"] = subdf.apply(lambda x: x.astype(str).sum(), axis=1)
subdf["hash"] = subdf["hash"].apply(lambda x: hashlib.shake_256(x.encode()).hexdigest(3))

subdf.to_csv(OUTPUT_FILE, sep="\t")
