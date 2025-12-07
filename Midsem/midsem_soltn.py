'''
NAME ASHIKA UTHAMAN
2311043
'''

import math
from My_Library import*

### Question 1 ###
lx,ly=LCG_rand(x0=0.1,a=1103515245, c=12345, m=32768,k=5,num_iter=10005)
# these Lx and Ly has numbers btw 0 and 1, since elipse a=2 and b=1, will replace Lx with random points btw 0 and 2, Ly remains same , hence we will solve in the 1st quadrant of ellipse part
lx_new=[]
for i in range(len(lx)):
    xnew=2*lx[i]
    lx_new.append(xnew)
# lx_new is matrix with random numbers btw 0 and 2

lin=[]#list containing values (x,y) with x**2 + 4*y**2<=4
lout=[]#list containing values (x,y) with x**2 + 4*y**2>4
area=[]#list of ellipse area obtained after each 20th iteration, area=pi*a*b=2*pi
iter_num=[]
for i in range(len(lx)):
    if (lx_new[i]**2 + (4*(ly[i]**2)))<=4:
        lin.append((lx_new[i],ly[i]))
    else:
        lout.append((lx_new[i],ly[i]))
    if (i+1)%20==0:
        area_upd=8*(len(lin)/(i+1)) ########## NOTE IN EXAM I ENTERED 4 TIMES THINKING ITS THE QUADRANT NUMBER, BUT IT IS THE AREA OF THE RWCTANGLE, MESSED UP....SO I GOT Pi AS SOLTN
        area.append(area_upd)
        iter_num.append(i+1)
Plot(iter_num, area, title='Area of ellipse vs number of iteration', xlabel='no.of iteration', ylabel='Ellipse area value', file_name='midsem_quest_1_plot.png')
print("Area value:" ,area_upd)


### Question two ###
def f(x):
    return (x-5)*math.exp(x) + 5
def df(x):
    return (x-5)*math.exp(x)+ math.exp(x)
x_solt=newton_raphson(f,df,3,10**(-4),10**(-4))
print("x_solt=the soltn obtained from newton raphson=",x_solt )

# lambda_m*T=hc/xk
h=6.626*(10**(-34))
c=3*(10**(8))
k=1.381*(10**(-23))
print("b=",(h*c)/(k*x_solt))


### Question3 ###
A1=read_matrix('mid_3_mat')
A2=read_matrix('mid_3_mat')
A3=read_matrix('mid_3_mat')
A4=read_matrix('mid_3_mat')
A5=read_matrix('mid_3_mat')



I=[[[1],[0],[0],[0],[0]],[[0],[1],[0],[0],[0]],[[0],[0],[1],[0],[0]],[[0],[0],[0],[1],[0]],[[0],[0],[0],[0],[1]]]
x1=LU_soltn(doolit_LU(A1),I[0])
x2=LU_soltn(doolit_LU(A2),I[1])
x3=LU_soltn(doolit_LU(A3),I[2])
x4=LU_soltn(doolit_LU(A4),I[3])
x5=LU_soltn(doolit_LU(A5),I[4])
Inverse=[]
Inverse.append(x1)
Inverse.append(x2)
Inverse.append(x3)
Inverse.append(x4)
Inverse.append(x5)
print(Inverse)# please mind that each entry in the matrix is inside a square bracket ########## NOTE...THIS IS ACTUALLY TRANSPOSE OF THE INVERSE MATRIX #########

### Question 4 ###
A=read_matrix('mid_4_mat')
b=read_matrix('mid_4_vect')
print("soltn is:", gauss_seidel(A,b))

################ RESULT ##################
'''
Area value: 3.1148
No. of iteration for Newton Raphson method is: 6
x_solt=the soltn obtained from newton raphson= 8.057602033435902e-14
b= 178637731069.28387
[[[-0.7079491478066469], [-0.19343159169227886], [0.021689019279128387], [0.27341203315637475], [0.7815264971593533]], [[2.5314333612741007], [0.31014249790444265], [0.3654652137468567], [-0.12992455993294275], [-2.875104777870915]], [[2.4311958647666922], [0.27946586569805326], [0.2861483654652136], [0.13161264785321758], [-2.678937319549222]], [[0.9665758591785414], [0.05771514389494271], [0.0505553227158424], [-0.14101355127130483], [-0.7011385861972617]], [[-3.9022771723945238], [-0.2941347676259662], [-0.2899203688181055], [0.44885675700847544], [4.233840923907981]]]
soltn is: [[1.4999998297596435], [-0.4999999999999992], [1.9999999999999996], [-2.499999914864037], [1.0000000000000004], [-0.9999999999957907]]

'''











