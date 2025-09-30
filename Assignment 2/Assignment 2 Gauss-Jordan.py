"""NAME: Ashika Uthaman
   Roll No.:2311043
"""

from My_Library import*
#let A be our Augmented matrix
A=[[0,2,5,1],[3,-1,2,-2],[1,-1,3,3]]
#V=[[1],[-2],[3]]

A1=read_matrix('asgn2_matA - Copy')
print(A1)
vect_A=read_matrix('asgn2_vect_A')
print(vect_A)
for i in range(len(vect_A)):
    A1[i].append(vect_A[i][0])
aug_A1=A1
print(aug_A1)


A2=read_matrix('asgn2_mat_n')
print(A2)
vect_A=read_matrix('asgn2_vect_n')
print(vect_A)
for i in range(len(vect_A)):
    A2[i].append(vect_A[i][0])
aug_A2=A2
print(aug_A2)

def Gauss_3(A):
    k=1
    if A[0][0]<A[1][0] or A[0][0]<A[2][0]: 
        if A[1][0]>A[2][0]:
            k=1
        else:
            k=2
    A[0],A[k]=A[k],A[0]
    r00=A[0][0]
    if A[0][0]!=1:
        for j in range(4):
            A[0][j]=A[0][j]/r00
    for j in range(1,3):
        if A[j][0]!=0:
            l1=A[0]
            l2=A[j]
            rj0=A[j][0]
            for k in range(4):
                A[j][k]=l2[k]-rj0*l1[k]
    #print(A)
    r11=A[1][1]
    if A[1][1]!=1:
        for j in range(4):
            A[1][j]=A[1][j]/r11
    for j in [0,2]:
        if A[j][1]!=0:
            l1=A[1]
            l2=A[j]
            r=A[j][1]
            for k in range(4):
                A[j][k]=l2[k]-r*l1[k]
    #print(A)
    r22=A[2][2]
    if A[2][2]!=1:
        for j in range(4):
            A[2][j]=A[2][j]/r22
    for j in [0,1]:
        if A[j][2]!=0:
            l1=A[2]
            l2=A[j]
            r=A[j][2]
            for k in range(4):
                A[j][k]=l2[k]-r*l1[k]
    print("The 3x3 row echelon matrix is",A)
    solt=[]
    for i in range(3):
        solt.append(A[i][3])
    return("The solutions are [x1,x2,x3]=",solt)
    
            
print(Gauss_3(aug_A1))    
        
#dim=len(A)
'''
def Gauss_Jordan(A,n=len(A)):
    #dimension of matrix without the augmented form is nxn
    #doing row wise pivoting and updating that row and setting the column elements in the pivoted row element as zero
    # only 1st row's element has to compared to other first elemnts and swap
    c=A[0][0]
    k=0
    for i in range(2,len(A)):
        r=A[1][0]
        if r<A[i][0]:
            r=A[i][0]
            k=i
    if A[0][0]<A[k][0]:
        A[0],A[k]=A[k],A[0]
    # row wise selected as pivoted row and doing elimination    
    for i in range(len(A)):
       
        piv_index=i #indicate the index of the pivot row
        
        rii=A[i][i]
        if A[i][i]!=1:
            for j in range(n+1):
                A[i][j]=A[i][j]/rii
        for k in range(n):
            if k!=piv_index:
                if A[k][i]!=0:
                    l1=A[i]
                    l2=A[k]
                    r=A[k][i]
                    for p in range(n+1):
                        A[k][p]=l2[p]-r*l1[p]
    print("The row echelon matrix is",A)
    solt=[]
    for i in range(n):
        solt.append(A[i][n])
    return("The solutions are [x1,x2,x3,...Xn]=",solt)
    
'''

print(Gauss_Jordan(aug_A2,n=len(aug_A2)))

############################################################    
#[[0.0, 2.0, 5.0], [3.0, -1.0, 2.0], [1.0, -1.0, 3.0]]
#[[1.0], [-2.0], [3.0]]
#[[0.0, 2.0, 5.0, 1.0], [3.0, -1.0, 2.0, -2.0], [1.0, -1.0, 3.0, 3.0]]
#[[1.0, -1.0, 4.0, 0.0, 2.0, 9.0], [0.0, 5.0, -2.0, 7.0, 8.0, 4.0], [1.0, 0.0, 5.0, 7.0, 3.0, -2.0], [6.0, -1.0, 2.0, 3.0, 0.0, 8.0], [-4.0, 2.0, 0.0, 5.0, -5.0, 3.0], [0.0, 7.0, -1.0, 5.0, 4.0, -2.0]]
#[[19.0], [2.0], [13.0], [-7.0], [-9.0], [2.0]]
#[[1.0, -1.0, 4.0, 0.0, 2.0, 9.0, 19.0], [0.0, 5.0, -2.0, 7.0, 8.0, 4.0, 2.0], [1.0, 0.0, 5.0, 7.0, 3.0, -2.0, 13.0], [6.0, -1.0, 2.0, 3.0, 0.0, 8.0, -7.0], [-4.0, 2.0, 0.0, 5.0, -5.0, 3.0, -9.0], [0.0, 7.0, -1.0, 5.0, 4.0, -2.0, 2.0]]
#The 3x3 row echelon matrix is [[1.0, 0.0, 0.0, -2.0], [0.0, 1.0, 0.0, -2.0], [0.0, 0.0, 1.0, 1.0]]
#('The solutions are [x1,x2,x3]=', [-2.0, -2.0, 1.0])
#The row echelon matrix is [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.7618170439978598], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.8962280338740136], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 4.051931404116158], [-0.0, -0.0, -0.0, 1.0, 0.0, 0.0, -1.6171308025395423], [-0.0, -0.0, -0.0, -0.0, 1.0, 0.0, 2.0419135385019125], [-0.0, -0.0, -0.0, -0.0, -0.0, 1.0, 0.1518324871559355]]
#('The solutions are [x1,x2,x3,...Xn]=', [-1.7618170439978598, 0.8962280338740136, 4.051931404116158, -1.6171308025395423, 2.0419135385019125, 0.1518324871559355])
                
            
            
            
            

    
     
      
            
                
                
                
        
                
        
        
        
        
