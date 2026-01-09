import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data=pd.read_csv("L2file.csv")

h=np.array(data.h)
L2=np.array(data.L2_error)

    
plt.title("L2-error versus h")
plt.plot(h,L2)
plt.xlabel("h")
plt.ylabel("L2 error")
plt.savefig("L2graph.jpg")



