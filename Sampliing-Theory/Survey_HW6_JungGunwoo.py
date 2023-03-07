import numpy as np
import pandas as pd
from sklearn import cluster

# Problem 24 in Section 6.9

pi = pd.DataFrame([0.31,0.2,0.14,0.03,0.01,0.31],
                  columns=["Pr"], index=["pi12", "pi13","pi14","pi23","pi24","pi34"])
'''
>>> pi
        Pr
pi12  0.31
pi13  0.20
pi14  0.14
pi23  0.03
pi24  0.01
pi34  0.31
'''
pi1 = pi["Pr"][0]+pi["Pr"][1]+pi["Pr"][2]
pi2 = pi["Pr"][0]+pi["Pr"][3]+pi["Pr"][4]
pi3 = pi["Pr"][1]+pi["Pr"][3]+pi["Pr"][5]
pi4 = pi["Pr"][2]+pi["Pr"][4]+pi["Pr"][5]
Pi = [pi1, pi2, pi3, pi4]
'''
>>> Pi
[0.65, 0.35, 0.54, 0.46]
'''

t =[2.5,2.0,1.1,0.5]

t_hat_HT = sum(t[i]/Pi[i] for i in range(0,len(t)))
# 12.684433119215727

table_ik = pd.DataFrame([(0.65, 0.35),(0.65, 0.54),(0.65, 0.46),(0.35, 0.54),(0.35, 0.46),(0.54, 0.46)],
                        columns=["pi_i", "pi_k"], index=["(1,2)","(1,3)","(1,4)","(2,3)","(2,4)","(3,4)",])
'''
>>> table_ik
       pi_i  pi_k
(1,2)  0.65  0.35
(1,3)  0.65  0.54
(1,4)  0.65  0.46
(2,3)  0.35  0.54
(2,4)  0.35  0.46
(3,4)  0.54  0.46
'''

table_ik["pi_ik"]=list(pi["Pr"])
table_ik["ti"]=[2.5,2.5,2.5,2.0,2.0,1.1]
table_ik["tk"]=[2.0,1.1,0.5,1.1,0.5,0.5]
table_ik["ti/pi_i"]=table_ik["ti"]/table_ik["pi_i"]
table_ik["tk/pi_k"]=table_ik["tk"]/table_ik["pi_k"]

'''
>>> table_ik
       pi_i  pi_k  pi_ik   ti   tk   ti/pi_i   tk/pi_k
(1,2)  0.65  0.35   0.31  2.5  2.0  3.846154  5.714286
(1,3)  0.65  0.54   0.20  2.5  1.1  3.846154  2.037037
(1,4)  0.65  0.46   0.14  2.5  0.5  3.846154  1.086957
(2,3)  0.35  0.54   0.03  2.0  1.1  5.714286  2.037037
(2,4)  0.35  0.46   0.01  2.0  0.5  5.714286  1.086957
(3,4)  0.54  0.46   0.31  1.1  0.5  2.037037  1.086957   
'''

V_hat_HT__t_hat_HT = sum((1-Pi[i])*(t[i]/Pi[i])**2 for i in range(0,len(Pi)))
+sum((1-table_ik["pi_i"][i]*table_ik["pi_k"][i]/table_ik["pi_ik"][i])*(table_ik["ti/pi_i"][i])*(table_ik["tk/pi_k"][i])
     for i in range(0,len(table_ik)))
# -130.9074063344015
'''
The estimator of variance of Horvitz-Thompson estimate is negative.
This indicates an "Unequal probability design"
'''

V_hat_SYG__t_hat_HT = 1/2 * sum(((table_ik["pi_i"][i]*table_ik["pi_k"][i]-table_ik["pi_ik"][i])/table_ik["pi_ik"][i])
    *(table_ik["ti"][i]/table_ik["pi_i"][i]-table_ik["tk"][i]/table_ik["pi_k"][i])**2 for i in range(0,len(table_ik)))
# 202.50028995281517



# Problem 26 in Section 6.9

cluster = pd.DataFrame([(1,5,"3,4,5,6,2",20),(2,4,"7,4,7,7",25),(3,8,"7,2,9,4,5,3,2,6",38),(4,5,"2,5,3,6,8",24),(5,3,"9,7,5",21),],
                       columns=["psi,i","Mi","Values,yij","ti"])
cluster["psi_i"]=cluster["Mi"]/sum(cluster["Mi"])
cluster["pi_i"]=cluster["psi_i"]*2
cluster["a_i"]=cluster["psi_i"]*(1-cluster["psi_i"])/(1-cluster["pi_i"])
'''
>>> cluster
   psi,i  Mi       Values,yij  ti  psi_i  pi_i       a_i
0      1   5        3,4,5,6,2  20   0.20  0.40  0.266667
1      2   4          7,4,7,7  25   0.16  0.32  0.197647
2      3   8  7,2,9,4,5,3,2,6  38   0.32  0.64  0.604444
3      4   5        2,5,3,6,8  24   0.20  0.40  0.266667
4      5   3            9,7,5  21   0.12  0.24  0.138947
'''

def pi_ij(i,j):
    if i==j:
        return 0
    else:
            pi_ij = cluster["a_i"][i-1]/sum(cluster["a_i"]) * cluster["psi_i"][j-1]/(1-cluster["psi_i"][i-1])
            + cluster["a_i"][j-1]/sum(cluster["a_i"]) * cluster["psi_i"][i-1]/(1-cluster["psi_i"][j-1])
    return pi_ij
pi_table = pd.DataFrame([1,2,3,4,5],columns=["i"],index=["1","2","3","4","5"])
pi_table["j=1"]=[pi_ij(1,1),pi_ij(1,2),pi_ij(1,3),pi_ij(1,4),pi_ij(1,5)]
pi_table["j=2"]=[pi_ij(2,1),pi_ij(2,2),pi_ij(2,3),pi_ij(2,4),pi_ij(2,5)]
pi_table["j=3"]=[pi_ij(3,1),pi_ij(3,2),pi_ij(3,3),pi_ij(3,4),pi_ij(3,5)]
pi_table["j=4"]=[pi_ij(4,1),pi_ij(4,2),pi_ij(4,3),pi_ij(4,4),pi_ij(4,5)]
pi_table["j=5"]=[pi_ij(5,1),pi_ij(5,2),pi_ij(5,3),pi_ij(5,4),pi_ij(5,5)]
'''
>>> pi_table
   i       j=1       j=2       j=3       j=4       j=5
1  1  0.000000  0.068091  0.192926  0.090434  0.048549
2  2  0.068091  0.000000  0.147531  0.068091  0.036286
3  3  0.192926  0.147531  0.000000  0.192926  0.106617
4  4  0.090434  0.068091  0.192926  0.000000  0.048549
5  5  0.048549  0.036286  0.106617  0.048549  0.000000
'''
V_t_hat_HT = 0
for i in range(0,4):
    for j in range(i+1,5):
        V_t_hat_HT += 1/2 * (cluster["pi_i"][i]*cluster["pi_i"][j]-pi_table.iloc[i,:][j+1])\
        *((cluster["ti"][i]/cluster["pi_i"][i]-cluster["ti"][j]/cluster["pi_i"][j])**2)
# 121.53316028931401