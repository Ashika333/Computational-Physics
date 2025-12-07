"""
Name:Ashika Uthaman
Roll. No: 2311043
"""

from My_Library import*
import math

def f1(x1,x2,x3): 
    return x1**2+x2-37
def f2(x1,x2,x3):
    return x1-x2**2-5
def f3(x1,x2,x3):
    return x1+x2+x3-3

def df(x1,x2,x3):
    
    df1x1=2*x1
    df1x2=1
    df1x3=0
    df2x1=1
    df2x2=-2*x2
    df2x3=0
    df3x1=1
    df3x2=1
    df3x3=1
    l=[[df1x1,df1x2,df1x3],[df2x1,df2x2,df2x3],[df3x1,df3x2,df3x3]]
    return l

    
    
def g1(x1,x2,x3):#for finding x1
    return math.sqrt(37-x2)
def g2(x1,x2,x3):#for finding x2
    return math.sqrt(x1-5)
def g3(x1,x2,x3):#for finding x3
    return 3-x1-x2


print(fixed_point_multi(g1,g2,g3,x01=5,x02=5,x03=5,eps=10**(-6)))#6,2,1

print("The solutions obtained from Neton Raphson is:",newton_raphson_multi(f1,f2,f3,df,x01=5,x02=5,x03=5,eps=10**(-6)))
