"""
   Name: Ashika Uthaman
   Roll.No.:2311043
   Assignment 1a.
"""


import matplotlib.pyplot as plt 
#define x and y , as a python list 
def Plot(x, y , title='Sample Plot', xlabel='X-axis Label', ylabel='Y-axis Label',file_name='sample_plot.png'):
    plt.plot(x, y, marker='o', linestyle='', color='b',label='')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)  
    plt.savefig(file_name)
    plt.show()

## Random generator using pRNG ##
def random_gen(c,x0=0.1,k=5,num_iter=1000):
    l=[x0]
    x=x0
    for i in range(1,num_iter):
        x=c*x*(1-x)
        l.append(x)
    print(l)
    
    L2=[]
    L3=[]
    for i in range(len(l)):
        if i+k<=len(l)-1:
            L2.append(l[i])
            L3.append(l[i+k])
    print(l)
    print(L2)
    print(L3)
   
    return( Plot(L2, L3 , title=(c,k), xlabel='x_i', ylabel='X_(i+k)',file_name='sample_plot.png'))
     
#print(random_gen(2.98,x0=0.1,k=6,num_iter=1000))#correlated
#print(random_gen(3.9135,x0=0.1,k=7,num_iter=1000))#little uncorrelated
#print(random_gen(3.28,x0=0.1,k=6,num_iter=1000))#correlated
#print(random_gen(3.68,x0=0.1,k=9,num_iter=1000))# correlated
#print(random_gen(3.86,x0=0.1,k=10,num_iter=1000))# uncorrelated


## LCG Random Generator ##
def LCG_rand(x0,a=1103515245, c=12345, m=32768,k=5,num_iter=1000):
    l=[x0]
    x=x0
    for i in range(1,num_iter):
        x=(a*x+c)%m
        x_i=x/m
        l.append(x_i)
    print(l)
    
    L2=[]
    L3=[]
    for i in range(len(l)):
        if i+k<=len(l)-1:
            L2.append(l[i])
            L3.append(l[i+k])
    #print(l)
    #print(L2)
    #print(L3)
    Plot(L2, L3 , title='Sample Plot', xlabel='x_i', ylabel='X_(i+k)',file_name='sample_plot.png')
    return(L2,L3)
print(LCG_rand(x0=0.1,a=1103515245, c=12345, m=32768,k=5,num_iter=1000))
     

#############################################   
    
#print(random_gen(2.98,x0=0.1,k=6,num_iter=1000))
#print(random_gen(3.9135,x0=0.1,k=7,num_iter=1000))
#print(random_gen(3.28,x0=0.1,k=6,num_iter=1000))
#print(random_gen(3.68,x0=0.1,k=9,num_iter=1000))
#print(random_gen(3.86,x0=0.1,k=10,num_iter=1000))
#print(LCG_rand(x0=0.1,a=1103515245, c=12345, m=32768,k=5,num_iter=1000))

###############################################################

        
    
