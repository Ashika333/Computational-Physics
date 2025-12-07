"""
NAME: Ashika Uthaman
Roll No.2311043
"""

import math
from My_Library import*

### Question 1###
def f(x):
    return math.log(x/2)-math.sin((5*x)/2)

#i=0

print("The root found for the fnctn f(x) in the interval [1.5,3.0] by bisection is:", root_bisect(f,1.5,3.0))

print("The root for the function f(x) in the interval [1.5,3.0] by  Regula Falsi method is:",root_regula_falsi(f,1.5,3.0))


### Question 2 ###

def g(x):
    return -x-math.cos(x)
def bracketing(g,a,b):# it returns an interval where the root of the function exist in this interval
    #g implies g(x) which is the given fnctn
    beeta=0.3
    i=0
    if g(a)*g(b)<0:
        print("value of g(a)*g(b) after bracketing is:",g(a)*g(b))
        return([a,b])
    else:
        if abs(g(a))<abs(g(b)):
            a=a-(beeta*(b-a))
            return(bracketing(g,a,b))
        elif abs(g(a))>abs(g(b)):
            b=b+(beeta*(b-a))
            return(bracketing(g,a,b))
print("Bracketing for the function g(x):",bracketing(g,2,4))


###############################################################
### RESULT ###
'''
The root found for the fnctn f(x) in the interval [1.5,3.0] by bisection is: 2.6231406927108765
The root for the function f(x) in the interval [1.5,3.0] by  Regula Falsi method is: 2.6231403358555436
value of g(a)*g(b) after bracketing is: -6.201243156109751
Bracketing for the function g(x): [-1.7122000000000002, 4]


'''
###############################################################
                
    
            
   
                
             
            
    
        
    
    
                
            
    
    
        
    
    
    
