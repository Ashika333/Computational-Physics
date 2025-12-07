'''
Name: Ashika Uthaman
Roll.No. 2311043
'''

import math
from My_Library import*

### Question 1 ###
def h(x):
    return 3*x+math.sin(x)-math.exp(x)
def dh(x):#dh is derivative of fnctn h
    return 3+math.cos(x)-math.exp(x)
'''
From the result it can be understood that Newton Raphson converges much faster than the bisectio and the Regula falsi method. And regula Falsi converges faster than bisection method.
These can be understood from the iteration
'''

print("The solution by Newton Raphson is:",newton_raphson(h,dh,1,10**(-6),10**(-6)))
print("The root value by bisection:",root_bisect(h,-1.5,1.5))
print("The root value by Regula falsi is:",root_regula_falsi(h,-1.5,1.5))

### Question 2 ###
def h1(x):
    return(x**2 -(2*x)-3)
def gh1(x):
    # g(x) to be defined for this h1(x) is represented by gh1(x)
    return (((x**2)-3)/2)

print("The root from Fixed Point method:",fixed_point(gh1,1,10**(-6))) # in this problem i ve noted that for my g(x)=(x**2-3)/2, at x1=1 choice only im gwtting, for other tried values like 0,0.5,2,..im not getting any poutput, wh
    
#########################################################################################################
### RESULT ###
'''
No. of iteration for Newton Raphson method is: 5
The solution by Newton Raphson is: 0.36042170296019965
The no.of iteration from bisection: 22
The root value by bisection: 0.3604220151901245
No.of iteration for the regula falsi method is: 11
The root value by Regula falsi is: 0.36042180364913606
The root from Fixed Point method: -1.0
'''
###########################################################################################################


