'''
NAME: ASHIKA UTHAMAN
ROLL.NO. 2311043
ASSIGNMENT 9: LAGUERRE'S METHOD
'''
import math
from My_Library import*
'''
def polynomial(x,n,l):
    #n=degree of polynomial
    #l=list containing the co-efficients of each term in polynomial whoes power ranges from o to n
    # if poly is 3*x^0 + 7*x + 2*x^3 ==> l=[3,7,0,2]
    #n degree poly means n +1 terms
    poly=0
    for i in range(n+1):
        poly+=l[i]*(x**i)
    return poly

def first_dpoly(x,n,l):# To find the first derivative of our polynomial
    d1_poly=0
    for i in range(1,n+1):
        d1_poly+=l[i]*i*(x**(i-1))
    return d1_poly

def second_dpoly(x,n,l): #To find the second derivative
    d2_poly=0
    for i in range(2,n+1):
        d2_poly+=l[i]*i*(i-1)*(x**(i-2))
    return d2_poly
'''
    
'''
list1=[]
def laguerre(n,l,beta0,eps,delta): #try eps=10**-6 and delta=10**-4, len(l)=n+1
    #for storing the factors and should have n factors equal to degree of polynomial
    k=10
    x0=beta0
    if polynomial(x0,n,l)<=delta:
        list1.append(x0)
        alpha=x0
    else:
        
        while k>eps and polynomial(x0,n,l)>delta:
            G=(first_dpoly(x0,n,l))/(polynomial(x0,n,l))
            H=G**2-((second_dpoly(x0,n,l))/(polynomial(x0,n,l)))
            val1=G+math.sqrt((n-1)*(n*H-(G**2)))
            val2=G-math.sqrt((n-1)*(n*H-(G**2)))
            denom=val1
            if abs(val2)>abs(val1):
                denom=val2
            a=n/denom
            xstore=x0
            x0=x0-a
            k=abs(x0-xstore)
        list1.append(x0)
        alpha=x0

    if len(list1)!=n:
        l_coeff=[ l[i] for i in range(n,-1,-1)]#for inversing the list l
        l_store=[l_coeff[0]]
        c=0
        s=0
        for i in range(len(l)-1):#
            c=alpha*l_store[-1]
            s=c + l_coeff[i+1]
            l_store.append(s)
        m=len(l_store)
        l_next=l_store[0:m-1]# the coeffients are from the higher power to lower power
        l_new=[l_next[i] for i in range(m-2,-1,-1)]#list with coefficients from 0th power of x till the new degree
        #therefore the new polynomial will be of degree n-1
        #degree of the new polynomial obtaining is one less than the length of new list
        
        deg_new=len(l_new)-1
        return (laguerre(deg_new,l_new,beta0,eps,delta)) #this return value is the alpha
    if len(list1)==n:
        return("The list of factors are:",list1)
'''
'''
def laguerre(n,l,beta0,eps,delta): #try eps=10**-6 and delta=10**-4, len(l)=n+1
    #for storing the factors and should have n factors equal to degree of polynomial
    k=10
    x0=beta0
    list=[]
    deg=n
    lp=l#lp=list of powers of the polynomial to be solved to find the root
    for i in range(n):
        if polynomial(x0,deg,lp)<=delta:
            list.append(x0)
            alpha=x0
        else:
            while k>eps and polynomial(x0,deg,lp)>delta:
                G=(first_dpoly(x0,deg,lp))/(polynomial(x0,deg,lp))
                H=G**2-((second_dpoly(x0,deg,lp))/(polynomial(x0,deg,lp)))
                val1=G+math.sqrt((n-1)*(n*H-(G**2)))
                val2=G-math.sqrt((n-1)*(n*H-(G**2)))
                denom=val1
                if abs(val2)>abs(val1):
                    denom=val2
                a=n/denom
                xstore=x0
                x0=x0-a
                k=abs(x0-xstore)
            list.append(x0)
            alpha=x0

        l_coeff=[ lp[i] for i in range(deg,-1,-1)]#for inversing the list l
        l_store=[l_coeff[0]]
        c=0
        s=0
        for i in range(len(lp)-1):#
            c=alpha*l_store[-1]
            s=c + l_coeff[i+1]
            l_store.append(s)
        m=len(l_store)
        l_next=l_store[0:m-1]# the coeffients are from the higher power to lower power
        l_new=[l_next[i] for i in range(m-2,-1,-1)]#list with coefficients from 0th power of x till the new degree
        #therefore the new polynomial will be of degree n-1
        #degree of the new polynomial obtaining is one less than the length of new list

        deg=deg-1
        lp=l_new
        x0=beta0

    return list
'''        

    
            
print("For the 1st polynomial",laguerre(4,[6,1,-7,-1,1],1.3,10**(-6),10**(-5)))  ###***** IF X0=1.3, THE SOLTN IS PRINTING, 1.3,1.3,1.3,3 INITIAL VALUE IS HIGHLY INFLUENCING*****###
print("For the 2nd polynomial",laguerre(4,[4,0,-5,0,1],2.4,10**(-6),10**(-5)))
'''
for the first two polynomials in the question it is printing error as "    a=n/denom
ZeroDivisionError: float division by zero". I wasnt able to figure out the error within the allowed time
'''
print("For the 3rd polynomial",laguerre(5,[-4.5,13.5,0.5,-19.5,0,2],6,10**(-6),10**(-5)))

'''
Since the code is not working, I didnt add it ti library
'''
#####################################################
'''
For the 3rd polynomial ('The list of factors are:', [3.0, 0.5006925035068446, 0.49930703946676447])
'''
#####################################################    
    
    
    
        
            







             
