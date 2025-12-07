'''
Name: Ashika Uthaman
Roll.No. 2311043
'''


import math
from My_Library import*


### Question 1 ###


A1=[[4,-1,0,-1,0,0],[-1,4,-1,0,-1,0],[0,-1,4,0,0,-1],[-1,0,0,4,-1,0],[0,-1,0,-1,4,-1],[0,0,-1,0,-1,4]]
b1=[[2],[1],[2],[2],[1],[2]]

'''
def diag_dom(A):
    n=len(A)
    for i in range(n):
        if i+1<=(n-1):
            c=A[i][i]
            d=A[i+1][i]
            s=i+1
            for j in range(i+2,n):
                if A[j][i]>d:
                    d=A[j][i]
                    s=j
        if d>c:
            A[i],A[s]=A[s],A[i]
    return(A)
print(diag_dom(A2))
'''
#my cholesky code not working and i am not able to figure it out why, so i didnt add it in library
'''
def cholesky(A):
    n=len(A)
    L=[]
    for p in range(n):
        L.append([])
        
    print(L)
    
    for i in range(n):
        for j in range(n):
            if j==i:
                sum=0
                for k in range(0,i):
                    sum+=(L[i][k])**2
                L[i].append(math.sqrt(A[i][i]-sum))

            if j<i:
                L[i].append(0)
            if j>i:
                sum=0
                for k in range(0,i):
                    sum+=L[i][k]*L[k][j]
                L[i].append((1/L[i][i])*(A[i][j]-sum))
            
    return(L)
'''
#checking symmetry
def sym_matrix(A):
    n=len(A)
    k=1
    for i in range(n):
        for j in range(n):
            if j!=i:
                c=A[i][j]
                if A[j][i]!=c:
                    #
                    k=0
                    break
    if k==0:
        return("The matrix is not symmetric")
    if k==1:
        return("This matrix is symmetric matrix")
    
print(sym_matrix(A1))
print("The solution from gauss seidel method:",gauss_seidel(A1,b1))
print("soltn from cholesky",cholesky_soltn(cholesky_decomp(A1),b1))

### Question 2 ###
A2=[[11,3,0,1,2],[0,4,2,0,1],[3,2,7,1,0],[4,0,4,10,1],[2,5,1,3,13]]
b2=[[51],[15],[15],[20],[92],]
print("The solution from gauss seidel method:",gauss_seidel(A2,b2))
print("The solution from jacobi method:",jacobi(A2,b2))

######################################################
'''
This matrix is symmetric matrix
The solution from gauss seidel method: ('The solutin matrix is:', [[0.9999997530614102], [0.9999997892247294], [0.9999999100460266], [0.9999998509593769], [0.9999998727858708], [0.9999999457079743]])
The solution from gauss seidel method: ('The solutin matrix is:', [[2.979165086347139], [2.2155996761867414], [0.21128402698819163], [0.15231700827754785], [5.715033568811629]])
The D matrix from jacobi method: [[11, 0, 0, 0, 0], [0, 4, 0, 0, 0], [0, 0, 7, 0, 0], [0, 0, 0, 10, 0], [0, 0, 0, 0, 13]]
The L matrix from jacobi method: [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [3, 2, 0, 0, 0], [4, 0, 4, 0, 0], [2, 5, 1, 3, 0]]
The U matrix from jacobi method: [[0, 3, 0, 1, 2], [0, 0, 2, 0, 1], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 0]]
The solution from jacobi method: [[2.979165056795253], [2.2155993914854326], [0.2112838649649856], [0.15231675099384498], [5.715033407349234]]
'''
########################################################
            
            
     
     
            
