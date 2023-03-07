import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics

# Q4
Data = pd.DataFrame([(15.0,14.0), (9.5, 9.0), (14.2, 12.5), (20.5, 22.0), (6.7, 6.3), (9.8, 8.4), (25.7, 28.5), (12.6, 10.0), (15.1, 14.4), (30.9, 28.2), (7.3, 15.5), (28.6, 26.3), (14.7, 13.1), (20.5, 19.5), (10.9, 9.8)], columns=["x","y"], index=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15",])
'''
       x     y
1   15.0  14.0
2    9.5   9.0
3   14.2  12.5
4   20.5  22.0
5    6.7   6.3
6    9.8   8.4
7   25.7  28.5
8   12.6  10.0
9   15.1  14.4
10  30.9  28.2
11   7.3  15.5
12  28.6  26.3
13  14.7  13.1
14  20.5  19.5
15  10.9   9.8
'''
# 4-a
plt.scatter(Data["x"], Data["y"])
plt.show()

r_sq = statistics.correlation(Data["x"], Data["y"]) # 0.932
r_sq = np.corrcoef(Data["x"], Data["y"])[0][1]

Cov_xy = statistics.covariance(Data["x"], Data["y"])
Cov_xy = np.cov(Data["x"], Data["y"])[0][1]
'''
array([[57.3052381 , 52.04309524],
       [52.04309524, 54.41238095]])
'''
s2_x = statistics.variance(Data["x"]) # 57.3052381
s2_y = statistics.variance(Data["y"]) # 54.41238095

# 4-d
Beta1_hat = Cov_xy/s2_x # 0.9081
Beta1_hat = r_sq*np.sqrt(s2_y)/np.sqrt(s2_x) # 0.9081
Beta0_hat = np.mean(Data["y"]) - Beta1_hat*np.mean(Data["x"]) # 1.1814
t_reg = Beta0_hat + Beta1_hat*45000 # 40868.9862

n = 15
N = 2100
V_t_reg = (1-n/N)*(1/n)*s2_y*(1-r_sq**2)*((45000)/(45000/2100))**2
# 2086566.6139

# 4-e
# (1)
t_y = np.mean(Data["y"]) * N # 33250

Var_t_y = (N*np.sqrt((s2_y/n)*(1-n/N)))**2 # 15882974.0

# (2)


# 4-f
Data["d"] = Data["y"]-Data["x"]
'''
       x     y    d
1   15.0  14.0 -1.0
2    9.5   9.0 -0.5
3   14.2  12.5 -1.7
4   20.5  22.0  1.5
5    6.7   6.3 -0.4
6    9.8   8.4 -1.4
7   25.7  28.5  2.8
8   12.6  10.0 -2.6
9   15.1  14.4 -0.7
10  30.9  28.2 -2.7
11   7.3  15.5  8.2
12  28.6  26.3 -2.3
13  14.7  13.1 -1.6
14  20.5  19.5 -1.0
15  10.9   9.8 -1.1
'''
t_reg = Beta0_hat + Beta1_hat*np.sum(Data["x"]) # 243.1814

t_diff = N * (np.mean(Data["y"])+45000/2100-np.mean(Data["x"])) # 44370

V_t_diff = N**2 * (1-n/N)*(1/n)*sum((Data["d"][i]-np.mean(Data["d"]))**2 for i in range(0,15))/(n-1) # 2227614





# Q5
Data = pd.DataFrame([(2, 2100, 4320), (6, 1910, 4160), (7, 3200, 5790)], columns=["Division","Number of employees", "Total number of sick-leave days"])
'''
   Division  Number of employees  Total number of sick-leave days
0         2                 2100                             4320     
1         6                 1910                             4160     
2         7                 3200                             5790  
'''
n=3
N=8
# 5-a
Data["Average numver of sick leaves"] = Data["Total number of sick-leave days"]/Data["Number of employees"]
'''
Average numver of sick leaves
            2.057143
            2.178010
            1.809375
'''
y_bar_hat_psi = np.mean(Data["Average numver of sick leaves"]) # 2.0148

# 5-b
M0 = np.sum(Data["Number of employees"]) # 7210
V_y_bar_hat_psi = np.sqrt(1/(3*2) * np.sum((Data["Average numver of sick leaves"][i]-y_bar_hat_psi)**2 for i in range(0,3))) # 0.1085

# 5-d
# a
y_bar_hat = np.sum(Data["Total number of sick-leave days"])/np.sum(Data["Number of employees"]) # 1.9792
# b
M0 = 12950
V_y_bar_hat = np.sqrt((1-n/N)/(3*(M0/N))*np.sum((Data["Number of employees"][i]**2)*((Data["Average numver of sick leaves"][i]-y_bar_hat)**2) for i in range(0,3))/(n-1))
# 5.4778
