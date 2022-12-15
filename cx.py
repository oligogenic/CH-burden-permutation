import pandas as pd
import sys 

f = sys.argv[1]

rf = pd.read_excel(f)

rf.to_csv(f.split(".")[0]+".csv",index= None, header = True, sep="\t")