from math import *
import numpy as np
import matplotlib.pyplot as plt
d=5.
r=100.
x = np.linspace(-3,3,1000)
lt=[ pow(10,i) for i in x]
def psi(t):
    return (r+1)/(r*d+1)-2*exp(-t)*r*(d-1)/((r*d)**2-1)

result=[]
resultapp=[]
for i in lt:
    result.append(psi(i))
    resultapp.append(1/d)


    
plt.figure(1)
plt.plot(lt,result)
plt.plot(lt,resultapp)
plt.xscale('log')
plt.ylim(0,1)
plt.show()
