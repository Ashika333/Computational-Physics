# NAME: ASHIKA UTHAMAN
#ROLL NO:2311043


from My_Library import*

### Question 1 ###
A1=[[1,2,4],[3,8,14],[2,6,13]]
print("The original matrix A1 is:",A1)
A1_LU=doolit_LU(A1)
print("A1 after LU decomposition:",A1_LU)
U=[]
L=[]
n=len(A1_LU)
for p in range(n):
    U.append([])
for k in range(n):
    L.append([])
for i in range(n):
    for j in range(n):
        if j>=i:
            U[i].append(A1_LU[i][j])
        else:
            U[i].append(0)
print("U matrix of A1 is:",U)
for i in range(n):
    for j in range(n):
        if j<i:
            L[i].append(A1_LU[i][j])
        elif j==i:
            L[i].append(1)
        else:
            L[i].append(0)
print("L matrix of A1:",L)
verify_matrix=matrix_mul(L,U)
print("On multiplying L and U matrix:",verify_matrix)
     
### Question 2 ###
'''
when i read the file containing matrixes in the 2nd question(both A and b), I am not able to get the matrix seperately. I meant, even if I assigned A2 and A3 in two different text files,
whatever changes I am making in A2 is getting updated in A3, which I doesn't want. So I ended up in typing the matrix A in the 2nd qustn in the functional code below this commending part.
A2=read_matrix('asgn2_mat_n')
b2=read_matrix('asgn2_vect_n')
A3=read_matrix('asgn2_mat_n2')
b3=read_matrix('asgn2_vect_n')
'''
#The reason why I had to type these matrixes instead of reading them from a text file is mentioned in the above commending part
A2=[[1,-1,4,0,2,9],[0,5,-2,7,8,4],[1,0,5,7,3,-2],[6,-1,2,3,0,8],[-4,2,0,5,-5,3],[0,7,-1,5,4,-2]]
b2=[[19],[2],[13],[-7],[-9],[2]]
A3=[[1,-1,4,0,2,9],[0,5,-2,7,8,4],[1,0,5,7,3,-2],[6,-1,2,3,0,8],[-4,2,0,5,-5,3],[0,7,-1,5,4,-2]]

    
    
print("The original matrix A2 is:",A2)
#print("A2 after LU decomposition:",doolit_LU(A2))
#print("Now A2(after passing it once through the doolit_LU function) is:",A2)# A2 has changed to its LU decomposed form

#########TRYING NEW CODE FOR FORWD AND BCKWD ####################
def LU_soltn(A_lu,b):
    #A_lu the nxn matrix is the LU decomposed matrix given and b is also given,nx1 matrix
    n=len(A_lu)
    y=[b[0]]
    for i in range(1,n):
        sm=0
        for j in range(i):
            sm+=A_lu[i][j]*y[j][0]
        y_i=b[i][0]-sm
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
print("soltn is:",LU_soltn(doolit_LU(A2),b2))
    
       
    
################### MY RESULTS ######################
'''
The original matrix A1 is: [[1, 2, 4], [3, 8, 14], [2, 6, 13]]
A1 after LU decomposition: [[1, 2, 4], [3.0, 2.0, 2.0], [2.0, 1.0, 3.0]]
U matrix of A1 is: [[1, 2, 4], [0, 2.0, 2.0], [0, 0, 3.0]]
L matrix of A1: [[1, 0, 0], [3.0, 1, 0], [2.0, 1.0, 1]]
On multiplying L and U matrix: [[1, 2.0, 4.0], [3.0, 8.0, 14.0], [2.0, 6.0, 13.0]]
The original matrix A2 is: [[1, -1, 4, 0, 2, 9], [0, 5, -2, 7, 8, 4], [1, 0, 5, 7, 3, -2], [6, -1, 2, 3, 0, 8], [-4, 2, 0, 5, -5, 3], [0, 7, -1, 5, 4, -2]]
soltn is: [[-1.761817043997862], [0.8962280338740133], [4.051931404116158], [-1.6171308025395421], [2.041913538501913], [0.15183248715593525]]

'''
########################################################         
             
            

    
           
            
    
    
        
    
