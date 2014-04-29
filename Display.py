import matplotlib.pyplot as plt  # plt is needed for ploting the output.
import os
import tkFileDialog
from os import listdir
from os.path import isfile, join

mypath="/home/thibault/ISEM_M1_cpp_simulation/simulated_data"

mypath=tkFileDialog.askdirectory(initialdir=mypath)
os.chdir(mypath)

onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
onlyfiles.sort()

color=["black","blue","red","green"]
for i,files in enumerate(onlyfiles):
    title=files.split("_")
    number=title[0][0:1]
    if title[0][0:2].isdigit():
        number=title[0][0:2]
    plt.figure(number)
    with open(files, 'r') as f:
        comment=f.readline()
        comment=comment.split(" ")
        tau=[]
        estimation=[]
        lower_bound=[]
        upperbound=[]
        for line in f:
            if not line[0:2]=="//" or line[0]=="#":
                list_line=line.split(" ")
                tau.append(float(list_line[0]))
                estimation.append(float(list_line[2]))
                lower_bound.append(float(list_line[3]))
                upperbound.append(float(list_line[4]))
    
    plt.plot(tau,estimation, color=color[i%4], label=str(title[1])+" and "+str(title[2][:-4]), linewidth=2)
    plt.plot(tau,lower_bound, color=color[i%4],linewidth=0.9,linestyle="--")
    plt.plot(tau,upperbound, color=color[i%4],linewidth=0.9,linestyle="--")
    plt.title("The estimation of relatedness for "+str(number)+" infecting nematodes\n")
    plt.legend(loc="best")
    plt.xlim(0, 6)

plt.show()


