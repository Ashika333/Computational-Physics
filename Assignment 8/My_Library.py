import matplotlib.pyplot as plt
import math
## Code for finding sum of the elements in the list
def sum_list(l):
    s=0
    for i in range (len(l)):
        s+=l[i]
    return s

## Finding factorial of ##
def fact(n):
    if n==0:
        return 1
    else:
        return n*fact(n-1)

## for converting the file into matrix ##
def read_matrix(filename):
    with open(filename,'r') as f:
        matrix=[]
        for line in f:
            row=[float(num) for num in line.strip().split()]
            matrix.append(row)
    return matrix

## matrix multiplication ##
def matrix_mul(A,B):
    m=len(A) #No.of rows
    n=len(A[0]) #No.of columns
    #print("dimension of", A," is:",m,n)
    p=len(B) #No.of rows and p=n
    q=len(B[0]) #No.of columns
    #print("dimension of", B," is:",p,q)
    AB=[]
    for i in range(m):
        row=[]
        c=0
        for k in range(q):
            for j in range(n):
                c= c+ A[i][j]*B[j][k]
            row.append(c)
            c=0
        AB.append(row)
    return AB

#define x and y , as a python list 
def Plot(x, y , title='Sample Plot', xlabel='X-axis Label', ylabel='Y-axis Label',file_name='sample_plot.png'):
    plt.plot(x, y, marker='o', linestyle='', color='b',label='')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)  
    plt.savefig(file_name)
    plt.show()

## LCG Random Generator ##
def LCG_rand(x0,a=1103515245, c=12345, m=32768,k=5,num_iter=1000):
    l=[x0]
    x=x0
    for i in range(1,num_iter):
        x=(a*x+c)%m
        x_i=x/m
        l.append(x_i)
    #print(l)
    
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

#To make the Augmented form of matrix A and vector b
def make_augmented(A,b):
    for i in range(len(b)):
        A[i].append(b[i][0])
    return A

#Solving linear equation using Gauss Jordan
def Gauss_Jordan(A):
    #please enter A as the augmented matrix
    #dimension of matrix without the augmented form is nxn
    #doing row wise pivoting and updating that row and setting the column elements in the pivoted row element as zero
    # only 1st row's element has to compared to other first elemnts and swap
    n=len(A)
    c=A[0][0]
    k=0
    r=A[1][0]
    for i in range(2,len(A)):
       #
       if r<A[i][0]:
            r=A[i][0]
            k=i
    if A[0][0]<A[k][0]:
        A[0],A[k]=A[k],A[0]
    # row wise selected as pivoted row and doing elimination . But before that comparing it with elements in that column to get largest element and swapping accordingly    
    for i in range(len(A)):
        if i>0:
            for k in range(i+1,n):
                if A[i][i]<A[k][i]:
                    A[i],A[k]=A[k],A[i]
        #now we fixed A[i][i] such that it is the largest element in its entire column        
        piv_index=i #indicate the index of the pivot row
        
        rii=A[i][i]
        if rii==0:
            return ("Matrix is singular")
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
    #print("The row echelon matrix is",A)
    solt=[]
    for i in range(n):
        solt.append([A[i][n]])
    return(solt)

#solving linear equation using LU decomposition
def doolit_LU(A):
    n=len(A)
    
    for j in range(n):
        for i in range(n):
            
            if 1<=i<j+1: 
                su=0
                for k in range(0,i):
                    su+=A[i][k]*A[k][j]
                A[i][j]=A[i][j]-su
                
                
            if j+1<=i<n: 
                lu=0
                for k in range(0,j):
                    lu+=A[i][k]*A[k][j]
                A[i][j]=(A[i][j]-lu)/A[j][j]
               
    return(A)#now on calling this fnctn on A, A has rearranged to its LU decomp from

#continuation of doolit to find solution of linnear eqtn using fwd and bwd
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

#code for Jacobi
def jacobi(A,b):
    n=len(A)
    for i in range(n):
        if A[i][i]==0:
            return("can't solve by Jacobi")
    '''
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
    '''
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
    #print("No. of iteration is:",k)
    
    return(x)

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

#code for gauss_seidel
def gauss_seidel(A,b):
    n=len(A)
    x=[]
    for p in range(n):
        x.append([0])
    eps=10
    c=0
    while eps>=(10)**(-6):
        
        k=0
        for i in range(n):
            
            sum1=0
            sum2=0
            for j in range(0,i):
                sum1+=A[i][j]*x[j][0]
            for j in range(i+1,n):
                sum2+=A[i][j]*x[j][0]
            c=(1/A[i][i])*(b[i][0]-sum1-sum2)
            k+=(x[i][0]-c)**2
            x[i][0]=c
        eps=math.sqrt(k)
        #print(eps)

    return(x)

#Bisect method
i=0#since "global i" is in the fnctn inside library, when we call the fnctn in this libray it will search for "i" defined inside library file
def root_bisect(f,a,b):#such that a and b contains a single root within them & a<b
    global i
    eps=10**(-6)
    dell=10**(-6)
    c=0
    
    if abs(b-a)<eps:
        if f(a)<dell or f(b)<dell:
            print("The no.of iteration from bisection:",i)
            return ((a+b)/2)
    else:
        c=(a+b)/2
        if f(c)*f(a)<0:
            b=c
            i=i+1
            
            return (root_bisect(f,a,b))
                
        elif f(c)*f(b)<0:
            a=c
            i=i+1
            
            return (root_bisect(f,a,b))
#Regula Falsi
def root_regula_falsi(f,a,b):#considering that a single root lies btw a and b,where a<b
    eps=10**(-6)
    dell=10**(-6)
    k=10
    c_before=0
    i=0
    

    if abs(b-a)<eps:
        if f(a)<dell or f(b)<dell:
            print(i)
            return ((a+b)/2)
    else:
        while k>=eps:
            
            c_after=b-(((b-a)*f(b))/(f(b)-f(a)))
            if f(a)*f(c_after)<0:
                b=c_after
                k=abs(c_before-c_after)
                c_before=c_after
                i=i+1
            elif f(b)*f(c_after)<0:
                a=c_after
                k=abs(c_before-c_after)
                c_before=c_after
                i=i+1
        print("No.of iteration for the regula falsi method is:",i)
        return(c_after)

#Newton Raphson
def newton_raphson(f,df,x0,eps,dell):
    #f=input fnctn #df=derivative of functn f
    #x0=our initial guess   #eps= we set the value, generally used 10**(-6)
    # dell value equals eps
    k=10
    x_before=x0
    i=0
    while k>=eps or f(x_before)>=dell:
        x_after=x_before-(f(x_before)/df(x_before))
        k=abs(x_after-x_before)
        x_before=x_after
        i=i+1
        
    print("No. of iteration for Newton Raphson method is:",i)
    return(x_after)

#Fixed point method
def fixed_point(gf,x1,eps):
    #gf(x) is the x=g(x)function for the f(x)
    #x1=initial value being chosen
    #eps= we set the value, generally used 10**(-6)
    k=10
    x_before=x1
    while k>=eps:
        x_after=gf(x_before)
        k=abs(x_after-x_before)
        x_before=x_after
    return(x_after)

def fixed_point_multi(g1,g2,g3,x01=6,x02=2,x03=1,eps=10**(-6)):
    k=10
    l=0
    x_bef=[x01,x02,x03]
    n=len(x_bef)
    while k>eps:
        
        sk=0
        sk1=0
        g=[g1(x_bef[0],x_bef[1],x_bef[2]),g2(x_bef[0],x_bef[1],x_bef[2]),g3(x_bef[0],x_bef[1],x_bef[2])]
        for i in range (n):
            
            xi=x_bef[i]
            sk+=x_bef[i]**2
            x_bef[i]=g[i] #g1(x_bef[0],x_bef[1],x_bef[2])
            g[0]=g1(x_bef[0],x_bef[1],x_bef[2])
            g[1]=g2(x_bef[0],x_bef[1],x_bef[2])
            g[2]=g3(x_bef[0],x_bef[1],x_bef[2])
            sk1+=(x_bef[i]-xi)**2
        sk1=math.sqrt(sk1)
        sk=math.sqrt(sk)
        k=sk1/sk
        l=l+1
        print("The new x values are:",x_bef, l)
    print("No. of iteration:",l)
    return x_bef

'''
adding both below code which i couldnt cmplt on 10.09.2025. pls check that, can ignore the print statements commended within the function which helped me to correct the mistakes
'''
#Inverse finding
def inverse(A):
    n=len(A)
    I=[[1 if j==i else 0 for j in range(n)] for i in range(n)]
    #print("The identity using",I)
    p=0
    c=A[0][0]
    k=0
    r=A[1][0]
    for i in range(2,len(A)):
       #
       if r<A[i][0]:
            r=A[i][0]
            k=i
    if A[0][0]<A[k][0]:
        A[0],A[k]=A[k],A[0]
        p=1
    if p==1:
        I[0],I[k]=I[k],I[0]
    #print("hello row exchange",I)
    #print("same in A how looks",A)
    for i in range(len(A)):
        # row wise selected as pivoted row and doing elimination . But before that comparing it with elements in that column to get largest element and swapping accordingly    
        if i>0:
            for k in range(i+1,n):
                if A[i][i]<A[k][i]:
                    A[i],A[k]=A[k],A[i]
                    I[i],I[k]=I[k],I[i]
        #now we fixed A[i][i] such that it is the largest element in its entire column        
        piv_index=i #indicate the index of the pivot row
        
        rii=A[i][i]
        if rii==0:
            return ("Matrix is singular")
        
        rii=A[i][i]
        if A[i][i]!=1:
            for j in range(n):
                A[i][j]=A[i][j]/rii
                I[i][j]=I[i][j]/rii
        #print("pivot element in I :",I)
        #print("how in A",A)
        for k in range(n):
            if k!=piv_index:
                if A[k][i]!=0:
                    l1=A[i]
                    l1i=I[i]
                    l2=A[k]
                    l2i=I[k]
                    r=A[k][i]
                    #print("the r is:", r)
                    #print("test1",l1i)
                    #print("In A",l1)
                    #print("test2",l2i)
                    #print("in A",l2)
                    for p in range(n):
                        A[k][p]=l2[p]-r*l1[p]
                        I[k][p]=l2i[p]-r*l1i[p]
                    #print("hloo",I)
                    #print("how A",A)
        #print("hihi", I)
        #print("hihi A",A)
    return I
def newton_raphson_multi(f1,f2,f3,df,x01=5,x02=5,x03=5,eps=10**(-6)):
    m=0
    n=3
    k=10
    x_bef=[x01,x02,x03]
    
    f=[[f1(x01,x02,x03)],[f2(x01,x02,x03)],[f3(x01,x02,x03)]]
    while k>eps:
        J=[[0,0,0],[0,0,0],[0,0,0]]
        for i in range(n):
            for j in range(n):
                J[i][j]=(df(x_bef[0],x_bef[1],x_bef[2])[i][j])
        #print("J matrix",J)
        
        '''J_inv=inverse(J)'''### either we can use this step or gauss jordan
        #print("Inverse of J",J_inv)
        '''x_matr=matrix_mul(J_inv,f)''' ### or we can find x_matr without using an inverse fnctn, i.e, by using gauss jordan. x_matr=J-1.f ==> J.x_matr=f which
                                         ### can be solved by Gauss Jordan
        aug_J=make_augmented(J,f)
        x_matr=Gauss_Jordan(aug_J)
        s1=0
        s2=0
        for j in range(n):
            xj=x_bef[j]
            x_bef[j]=x_bef[j]-x_matr[j][0]
            s1+=x_bef[j]**2
            s2+=(xj-x_bef[j])**2
        f=[[f1(x_bef[0],x_bef[1],x_bef[2])],[f2(x_bef[0],x_bef[1],x_bef[2])],[f3(x_bef[0],x_bef[1],x_bef[2])]]
        s2=math.sqrt(s2)
        s1=math.sqrt(s1)
        k=s2/s1
        m=m+1
        #print("new x values",x_bef)
    print("The no. of iteration for newton raphson:",m)
    return (x_bef)
