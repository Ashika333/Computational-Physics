'''
NAME: ASHIKA UTHAMAN
ROLL NO. 2311043
ASSIGNMENT 4
'''
import math
from My_Library import*

# rem for jaccob for a particular iterartion time, store all the new updating x^k values in a seperate liste and once after all of them are updated and then only replace all the previous values with the new values, to avoid overwriting.
A0=[[4,12,-16],[12,37,-43],[-16,-43,98]]
A1=[[4,1,1,1],[1,3,-1,1],[1,-1,2,0],[1,1,0,2]]
b1=[[3],[3],[1],[3]]
A2=[[4,1,1,1],[1,3,-1,1],[1,-1,2,0],[1,1,0,2]]
b2=[[3],[3],[1],[3]]

#cholesky decompostion to make the decomposed matrix
def cholesky_decomp(A):#Aco-eff matrix for the given linear eqtn
    n=len(A)
    for i in range(n):
        for j in range(n):
            if j==i:
                sm=0
                for k in range(0,i):
                    sm+=A[i][k]**2
                sm=A[i][i]-sm
                A[i][j]=math.sqrt(sm)
                
            if j>i:
                sm=0
                for k in range(0,i):
                    sm+=A[i][k]*A[k][j]
                l_ij=A[i][j]-sm
                A[i][j]=l_ij/A[i][i]
                A[j][i]=A[i][j]
    return(A)
print("chole decomp matrix form",cholesky_decomp(A0))

def cholesky_soltn(A_lu,b):
    #A_lu the nxn matrix is the cholesky decomposed matrix given and b is also given,nx1 matrix
    n=len(A_lu)
    y_0=b[0][0]/A_lu[0][0]
    y=[[y_0]]
    for i in range(1,n):
        sm=0
        for j in range(i):
            sm+=A_lu[i][j]*y[j][0]
        y_i=(b[i][0]-sm)/A_lu[i][i]
        y.append([y_i])#y is an nx1 matrix form
    #print("y is",y)
    
    x=[[0] for i in range(n)]
    x[n-1][0]=y[n-1][0]/A_lu[n-1][n-1]
    for i in range(n-2,-1,-1):
        sm=0
        for j in range(i+1,n):
            sm+=A_lu[i][j]*x[j][0]
        x_i=(y[i][0]-sm)/A_lu[i][i]
        x[i][0]=x_i
    return x
print("cholesky soltn",cholesky_soltn(cholesky_decomp(A1),b1))

'''
def cholesky(A):
    n=len(A)
    L=[]
    for p in range(n):
        L.append([])
        
    for i in range(n):
        for j in range(n):
            if j==i:
                sum1=0
                for k in range(0,i):
                    sum1+=(L[i][k])**2
                L[i].append(math.sqrt(A[i][i]-sum1))
                
            if j>i:
                sum=0
                for k in range(0,i):
                    sum+=L[i][k]*L[k][j]
                L[i].append((1/L[i][i])*(A[i][j]-sum))
                L[j].append((1/L[i][i])*(A[i][j]-sum))
                #print('edo',(1/L[i][i])*(A[i][j]-sum))
            
    return(L)
'''
#print(cholesky(A0))
'''
def jacobi(A,b):
    n=len(A)
    for i in range(n):
        if A[i][i]==0:
            return("can't solve by Jacobi")
    
    D=[]
    for p in range(n):
        D.append([])
    
    for i in range(n):
        for j in range(n):
            if j==i:
                D[i].append(A[i][j])
            else:
                D[i].append(0)
    print(D)
    L=[]
    U=[]
    for p in range(n):
        L.append([])
    for i in range(n): 
        for j in range(n):
            if j<i:
                L[i].append(A[i][j])
            else:
                L[i].append(0)
    print(L)
    
    for k in range(n):
        U.append([])
    for i in range(n):
        for j in range(n):
            if j>i:
                U[i].append(A[i][j])
            else:
                U[i].append(0)
           
    print(U)
    #print(A)
    
    x=[]
    store=[]
    for p in range(n):
        x.append([0])
    print(x)
    eps=10
    k=0
    while eps>=10**(-6):
        for i in range(n):
            sum=0
            for j in range(n):
                if j!=i:
                    sum+=A[i][j]*x[j][0]
            store.append((1/A[i][i])*(b[i][0]-sum))
        sum1=0
        for i in range(n):
            sum1+=(x[i][0]-store[i])**2
        eps=math.sqrt(sum1)
        print(eps)
        for p in range(n):
            x[p][0]=store[p]
        store=[]#missed
        k+=1
        print(k)
    #print(k)
        
        
        
    return(x)

def jacobi(A,b):
    n=len(A)
    for i in range(n):
        if A[i][i]==0:
            return("can't solve by Jacobi")

    D=[]
    for p in range(n):
        D.append([])
    
    for i in range(n):
        for j in range(n):
            if j==i:
                D[i].append(A[i][j])
            else:
                D[i].append(0)
    print("The D matrix from jacobi method:",D)
    L=[]
    U=[]
    for p in range(n):
        L.append([])
    for i in range(n): 
        for j in range(n):
            if j<i:
                L[i].append(A[i][j])
            else:
                L[i].append(0)
    print("The L matrix from jacobi method:",L)
    
    for k in range(n):
        U.append([])
    for i in range(n):
        for j in range(n):
            if j>i:
                U[i].append(A[i][j])
            else:
                U[i].append(0)
           
    print("The U matrix from jacobi method:",U)
    #print(A)

    #x=[[0]for i in range(n)]
    x=[]
    store=[]
    for p in range(n):
        x.append([0])
    #print(x)
    eps=10
    k=0
    while eps>=10**(-6):
        for i in range(n):
            sum=0
            for j in range(n):
                if j!=i:
                    sum+=A[i][j]*x[j][0]
            store.append((1/A[i][i])*(b[i][0]-sum))
        sum1=0
        for i in range(n):
            sum1+=(x[i][0]-store[i])**2
        eps=math.sqrt(sum1)
        #print(eps)
        for p in range(n):
            x[p][0]=store[p]
        store=[]#missed
        k+=1
        #print(k)
    
        
        
    return(x)
'''
print(jacobi(A2,b2))


'''
[[0, 0, 0], [0, 0, 0], [0, 0, 0]]
[[2.0, 5.744562646538029, 7.810249675906654], [2.0, 5.744562646538029, 7.810249675906654], [2.0, 5.744562646538029, 7.810249675906654]]
[[4, 0, 0, 0], [0, 3, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]
[[0, 0, 0, 0], [1, 0, 0, 0], [1, -1, 0, 0], [1, 1, 0, 0]]
[[0, 1, 1, 1], [0, 0, -1, 1], [0, 0, 0, 0], [0, 0, 0, 0]]
[[0], [0], [0], [0]]
[[0.75], [1.0], [0.5], [1.5]]
'''
    
                    
                    

                
    
