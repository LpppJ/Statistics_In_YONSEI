from random import sample
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
from statistics import LinearRegression

### Problem 9 in Section 3.9 of the textbook
table = pd.read_csv("agstrat.csv",index_col='county')
table.head()
n=300
N=3078

plt.scatter(table["acres87"],table["farms92"])
plt.show()
plt.scatter(table["acres87"],table["largef92"])
plt.show()
plt.scatter(table["acres87"],table["smallf92"])
plt.show()
plt.scatter(table["farms92"],table["largef92"])
plt.show()
plt.scatter(table["farms92"],table["smallf92"])
plt.show()
plt.scatter(table["largef92"],table["smallf92"])
plt.show()

acres87 = table["acres87"]
farms92 = table["farms92"]
smallf92 = table["smallf92"]
largef92 = table["largef92"]
region = table["region"]

stratum = pd.DataFrame([(220,21),(1054,103),(1382,135),(422,41)],
                       columns=("Nh","nh"), index=("NE","NC","S","W"))
nh=stratum["nh"]
Nh=stratum["Nh"]
H=len(stratum) #4

# acres87
strNE_acres87 = []
strNC_acres87 = []
strS_acres87 = []
strW_acres87 = []
for i in range(0,n):
    if region[i]=="NE":
        strNE_acres87.append(acres87[i])
    if region[i]=="NC":
        strNC_acres87.append(acres87[i])
    if region[i]=="S":
        strS_acres87.append(acres87[i])
    if region[i]=="W":
        strW_acres87.append(acres87[i])
str_acres87=[strNE_acres87,strNC_acres87,strS_acres87,strW_acres87]
sum_str_acres87 = [np.sum(str_acres87[0]),np.sum(str_acres87[1]),np.sum(str_acres87[2]),np.sum(str_acres87[3])]
stratum["sum_str_acres87"]=sum_str_acres87
stratum["mean_acres87"]=stratum["sum_str_acres87"]/stratum["nh"]

Sh_square = [0,0,0,0]
for i in range(0,H):
    for j in range(0,len(str_acres87[i])):
        Sh_square[i]+=((str_acres87[i][j]-stratum["mean_acres87"][i])**2) / (nh[i]-1)
stratum["Sh_square"]=Sh_square

var = [0,0,0,0]
for i in range(0,H):
    var[i]=(1-nh[i]/Nh[i])*(nh[i]/N)**2*(Sh_square[i]/nh[i])
stratum['Varhat_ystrbar'] = var

Vhat_ystrbar = stratum["Varhat_ystrbar"].sum() # 2,516,296.2739
SE_acres87_str_bar = np.sqrt(Vhat_ystrbar) # 1586.2838

# farms92
strNE_farms92 = []
strNC_farms92 = []
strS_farms92 = []
strW_farms92 = []
for i in range(0,n):
    if region[i]=="NE":
        strNE_farms92.append(farms92[i])
    if region[i]=="NC":
        strNC_farms92.append(farms92[i])
    if region[i]=="S":
        strS_farms92.append(farms92[i])
    if region[i]=="W":
        strW_farms92.append(farms92[i])
str_farms92=[strNE_farms92,strNC_farms92,strS_farms92,strW_farms92]
sum_str_farms92 = [np.sum(str_farms92[0]),np.sum(str_farms92[1]),np.sum(str_farms92[2]),np.sum(str_farms92[3])]
stratum["sum_str_farms92"]=sum_str_farms92
stratum["mean_farms92"]=stratum["sum_str_farms92"]/stratum["nh"]

Sh_square = [0,0,0,0]
for i in range(0,H):
    for j in range(0,len(str_farms92[i])):
        Sh_square[i]+=((str_farms92[i][j]-stratum["mean_farms92"][i])**2) / (nh[i]-1)
stratum["Sh_square"]=Sh_square

var = [0,0,0,0]
for i in range(0,H):
    var[i]=(1-nh[i]/Nh[i])*(nh[i]/N)**2*(Sh_square[i]/nh[i])
stratum['Varhat_ystrbar'] = var

Vhat_ystrbar = stratum["Varhat_ystrbar"].sum() # 5.5997
SE_farms92_str_bar = np.sqrt(Vhat_ystrbar) # 2.3664

# largef92
strNE_largef92 = []
strNC_largef92 = []
strS_largef92 = []
strW_largef92 = []
for i in range(0,n):
    if region[i]=="NE":
        strNE_largef92.append(largef92[i])
    if region[i]=="NC":
        strNC_largef92.append(largef92[i])
    if region[i]=="S":
        strS_largef92.append(largef92[i])
    if region[i]=="W":
        strW_largef92.append(largef92[i])
str_largef92=[strNE_largef92,strNC_largef92,strS_largef92,strW_largef92]
sum_str_largef92 = [np.sum(str_largef92[0]),np.sum(str_largef92[1]),np.sum(str_largef92[2]),np.sum(str_largef92[3])]
stratum["sum_str_largef92"]=sum_str_largef92
stratum["mean_largef92"]=stratum["sum_str_largef92"]/stratum["nh"]

Sh_square = [0,0,0,0]
for i in range(0,H):
    for j in range(0,len(str_largef92[i])):
        Sh_square[i]+=((str_largef92[i][j]-stratum["mean_largef92"][i])**2) / (nh[i]-1)
stratum["Sh_square"]=Sh_square

var = [0,0,0,0]
for i in range(0,H):
    var[i]=(1-nh[i]/Nh[i])*(nh[i]/N)**2*(Sh_square[i]/nh[i])
stratum['Varhat_ystrbar'] = var

Vhat_ystrbar = stratum["Varhat_ystrbar"].sum() # 0.1203
SE_largef92_str_bar = np.sqrt(Vhat_ystrbar) # 0.3469

# smallf92
strNE_smallf92 = []
strNC_smallf92 = []
strS_smallf92 = []
strW_smallf92 = []
for i in range(0,n):
    if region[i]=="NE":
        strNE_smallf92.append(smallf92[i])
    if region[i]=="NC":
        strNC_smallf92.append(smallf92[i])
    if region[i]=="S":
        strS_smallf92.append(smallf92[i])
    if region[i]=="W":
        strW_smallf92.append(smallf92[i])
str_smallf92=[strNE_smallf92,strNC_smallf92,strS_smallf92,strW_smallf92]
sum_str_smallf92 = [np.sum(str_smallf92[0]),np.sum(str_smallf92[1]),np.sum(str_smallf92[2]),np.sum(str_smallf92[3])]
stratum["sum_str_smallf92"]=sum_str_smallf92
stratum["mean_smallf92"]=stratum["sum_str_smallf92"]/stratum["nh"]

Sh_square = [0,0,0,0]
for i in range(0,H):
    for j in range(0,len(str_smallf92[i])):
        Sh_square[i]+=((str_smallf92[i][j]-stratum["mean_smallf92"][i])**2) / (nh[i]-1)
stratum["Sh_square"]=Sh_square

var = [0,0,0,0]
for i in range(0,H):
    var[i]=(1-nh[i]/Nh[i])*(nh[i]/N)**2*(Sh_square[i]/nh[i])
stratum['Varhat_ystrbar'] = var

Vhat_ystrbar = stratum["Varhat_ystrbar"].sum() # 0.4904
SE_smallf92_str_bar = np.sqrt(Vhat_ystrbar) # 0.7003

'''
SRS : Chapter2 #15
SE_acres87_bar = 18,913.6662
SE_farms92_bar = 22.0625
SE_largef92_bar = 3.9904
SE_smallf92_bar = 3.6375

Str : Chapter3 #9
SE_acres87_str_bar = 1586.2838
SE_farms92_str_bar = 2.3664
SE_largef92_str_bar = 0.3469
SE_smallf92_str_bar = 0.7003

The SE of the 'Str estimator' is smaller than the SE of 'SRS estimator'
'''
stratum.to_csv("stratum.csv")




### Problem 2 in Section 4.8 of the textbook
pop = pd.DataFrame([(13,10),(7,7),(11,13),(12,17),(4,8),(3,1),(11,15),(3,7),(5,4)],
columns=("x","y"),
index=("1","2","3","4","5","6","7","8","9"))
'''
pop
    x   y
1  13  10
2   7   7
3  11  13
4  12  17
5   4   8
6   3   1
7  11  15
8   3   7
9   5   4
'''
# 2(a) : tx, ty, Sx, Sy, R, B
x = pop["x"]
y = pop["y"]
tx = sum(x) #69
ty = sum(y) #82
N=len(pop) #9

SumofSquareX = 0
for i in range(0,N):
    SumofSquareX += (x[i]-np.mean(x))**2
Sx=np.sqrt(SumofSquareX/(N-1))
Sx # 4.09

SumofSquareY = 0
for i in range(0,N):
    SumofSquareY += (y[i]-np.mean(y))**2
Sy=np.sqrt(SumofSquareY/(N-1))
Sy # 5.1828

SumofXY = 0
for i in range(0,N):
    SumofXY += (x[i]-np.mean(x))*(y[i]-np.mean(y))
R = SumofXY/((N-1)*Sx*Sy)
R # 0.8152

B = ty/tx
B # 1.1884

# 2(b)
size = ("1","2","3","4","5","6","7","8","9")
combination = (list(map(''.join, it.combinations(size,3))))
n=3
sample1 = list()
for i in range(0,84):
    sample1.append(combination[i][0])
sample2 = list()
for i in range(0,84):
    sample2.append(combination[i][1])
sample3 = list()
for i in range(0,84):
    sample3.append(combination[i][2])

index_list = list(range(1,len(combination)+1))

table = pd.DataFrame({'Samples':combination},index=index_list)
table["sample1"]=sample1
table["sample2"]=sample2
table["sample3"]=sample3

xbar=list()
for i in range(1,85):
    xbar.append( (x[table["sample1"][i]] + x[table["sample2"][i]] + x[table["sample3"][i]] )/3 )
table["xbar"]=xbar

ybar=list()
for i in range(1,85):
    ybar.append( (y[table["sample1"][i]] + y[table["sample2"][i]] + y[table["sample3"][i]] )/3 )
table["ybar"]=ybar

table["Bhat"]=table["ybar"]/table["xbar"]
table["tsrs"]=table["ybar"]*N
table["tyr"]=table["Bhat"]*tx

table.to_csv("table.csv")

# 2(c)
plt.hist(table["tyr"])
plt.show()
plt.hist(table["tsrs"])
plt.show()
'''
If we compare two histograms,
we observe that the histograms are similar.
We also observe that ratio estimator is a little less spread out.
'''

# 2(d)
np.mean(table["tyr"]) #82.7924
np.var(table["tyr"]) 
np.mean(table["tsrs"]) # 82
bias_of_tyr = abs(np.mean(table["tyr"])-np.mean(table["tsrs"])) # 0.7924

# 2(e)
bias_of_yrbar = (1-n/N)*1/n/np.mean(x)*(B*Sx**2-R*Sx*Sy) # 0.07578
'''
Thus, bias_of_tyr approximates N*bias_of_yrbar = 9*0.07578
'''


### Problem 11 in Section 4.8 of the textbook
table = pd.read_csv("counties.csv")
table.head()

N=3141
n=100

# 11(a)
plt.hist(table["physician"])
plt.show()

# 11(b)
x=totpop = table["totpop"]
y=physician = table["physician"]

Total_physicians = N*np.mean(physician) # 3,141 * 297.17 =933,410.97

SumofSquareY = 0
for i in range(0,100):
    SumofSquareY += (physician[i]-np.mean(physician))**2
Sy=np.sqrt(SumofSquareY/(n-1)) # std of Sample Physician = 1,591.8705

SE_Nybar = N*np.sqrt((1-n/N)*(Sy**2)/n) # standard error of Nybar = 491,983

# 11 (c)
plt.scatter(totpop,physician)
plt.show()
'''
Above plot say that Ratio estimation is the more appropriate for this data set
Ratio estimation works best if the data are well fit by a straight line through the origin.
However, since the data are straight, the reg estimator also can be considered
'''

# 11 (d)
# Ratio estimator
Bhat = np.mean(physician)/np.mean(totpop) # 0.002507
tx = 255077536
tyr = Bhat*tx # 639506.0240

s_e_square =0
for i in range(0,n):
    s_e_square += (physician[i]-Bhat*totpop[i])**2
s_e_square = s_e_square/(n-1) # 172267.8739
s_e = np.sqrt(s_e_square) # 415.05

SE_tyr = N*np.sqrt(1-n/N)*tx/np.mean(x)*s_e/np.sqrt(n)
SE_tyr # 276047617.8333

# Regression estimator
model = LinearRegression(x,y).fit(totpop.values.reshape(-1,1),physician)




