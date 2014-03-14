import matplotlib.pyplot as plt  # plt is needed for ploting the output.
import os
import tkFileDialog
from os import listdir
from os.path import isfile, join

mypath="C:\Users\Thibault\Desktop\ISEM_M1_Cpp_Simulations\simulated_data"

mypath=tkFileDialog.askdirectory(initialdir=mypath)
os.chdir(mypath)

onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]

for i,files in enumerate(onlyfiles):
    plt.figure(i)
    title=files.split("_")
    with open(files, 'r') as f:
        comment=f.readline()
        comment=comment.split(" ")
        pop=[]
        estimation=[]
        lower_bound=[]
        upperbound=[]
        relatedness=[]
        for line in f:
            if not line[0:2]=="//" or line[0]=="#":
                list_line=line.split(" ")
                pop.append(float(list_line[0]))
                estimation.append(float(list_line[1]))
                lower_bound.append(float(list_line[2]))
                upperbound.append(float(list_line[3]))
                relatedness.append(float(list_line[4]))
    
    plt.plot(pop,estimation, color="black", label=comment[1], linewidth=2)
    plt.plot(pop,lower_bound, color="blue", label=comment[2], linewidth=1,linestyle="--")
    plt.plot(pop,upperbound, color="blue", label=comment[3], linewidth=1,linestyle="--")
    plt.plot(pop,relatedness, color="red", label=comment[4], linewidth=1.5)
    plt.title("The estimation of conditional expectation \n for "+str(title[0])+", "+str(title[1])+" and "+str(title[2]))
    plt.legend(loc="best")

plt.show()
