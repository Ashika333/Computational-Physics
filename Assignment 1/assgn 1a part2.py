"""
   Name: Ashika Uthaman
   Roll.no.: 2311043
"""


import matplotlib.pyplot as plt
import math
import numpy as np
from My_Library import *


###   QUESTION 3: Determining Pi value  ###
Lx,Ly= LCG_rand(x0=0.1,a=1103515245, c=12345, m=32768,k=5,num_iter=10005) #Lx gives the list with all random x_i values and Ly is list with X_i+k


lin=[] #list containing values with x**2 + y**2 <=1
lout=[]#list containing values with x**2+y**2>1
Pi=[]# list containing pi value we got for each 20th repetative iteration performed
iter_num=[]
for i in range(len(Lx)):
    if (Lx[i]**2 + Ly[i]**2)<=1:
        lin.append((Lx[i],Ly[i]))
    else:
        lout.append((Lx[i],Ly[i]))
    if (i+1)%20==0:
        Pi_upd=4*(len(lin)/(i+1))
        Pi.append(Pi_upd)
        iter_num.append(i+1)
# the following plot is for Pi vs no. of iterations
Plot(iter_num, Pi, title='Pi vs number of iteration', xlabel='no.of iteration', ylabel='Pi value', file_name='assgn1_pi_plot.png')    
#print("last pi value is",Pi_upd)# this gives the last 10,000th pi value 

l1=Pi[199:499]     #from the plot btw pi vs no.of iteration i have selected the saturated portion to calculate the average
Pi_value=(sum_list(l1))/len(l1)
print("The average Pi value is:", Pi_value)


###   QUESTION 4   ###

l1,l2=LCG_rand(x0=0.1,a=1103515245, c=12345, m=32768,k=5,num_iter=5005)
loglist=[] #list containing -ln(x) values of the random numbers in the range [0,1]
for i in range(len(l1)):
    loglist.append(-np.log(l1[i]))
plt.hist(loglist, bins=40, color='skyblue', edgecolor='black')
plt.savefig('histogram.png')
plt.show()


####################################################################################################################
#The average Pi value is: 3.1091465784948613
#####################################################################################################################
    

