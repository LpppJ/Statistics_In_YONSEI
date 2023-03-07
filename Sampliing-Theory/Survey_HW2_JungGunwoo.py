from matplotlib.pyplot import axis
import numpy as np
import pandas as pd
from pyparsing import col

# 2
# Obtain the standard error of the estimate of the population total of the elephants heavier than 8 tons
 
# (a) Under SRS, n=30
data = pd.DataFrame([(80,0.8),(120,0.4)], columns=["Nh","phU"], index=["Male","Female"])
'''
         Nh  phU
Male     80  0.8
Female  120  0.4
'''
# p = (# of elephant population > 8t) / (# elephant population)
p = (data["Nh"][0]*data["phU"][0]+data["Nh"][1]*data["phU"][1]) / (data["Nh"][0]+data["Nh"][1])
# p = 0.56

N = 200
n = 30

V_t_SRS = N**2*(N-n) / (n*(N-1)) * p*(1-p) # 280.657
SE_t_SRS = np.sqrt(V_t_SRS) # 16,753

# (b) Under Stratified SRS
n_male = n * data["Nh"][0]/(data["Nh"][0]+data["Nh"][1])
n_female = n * data["Nh"][1]/(data["Nh"][0]+data["Nh"][1])
n_str = [n_male, n_female] # 12, 18

var_male = data["Nh"][0]/(data["Nh"][0]-1) * data["phU"][0]*(1-data["phU"][0]) # 0.1620
var_female = data["Nh"][1]/(data["Nh"][1]-1) * data["phU"][1]*(1-data["phU"][1]) # 0.2420
var_str = [var_male, var_female]

V_t_str = sum((1-n_str[i]/data["Nh"][i]) * data["Nh"][i]**2 * var_str[i]/n_str[i] for i in range(0,2)) # 238.0229
SE_t_str = np.sqrt(V_t_str) # 15.428



# 4. Problem 3.6
data = pd.DataFrame([(35000,2,504),(45000,1,324),(10000,1,72)], columns=["Nh","Rh","nh"], index=["Houses","Apartments","Condos"])
data["NhRh"] = data["Nh"]*data["Rh"]
'''
               Nh  Rh   nh   NhRh
Houses      35000   2  504  70000
Apartments  45000   1  324  45000
Condos      10000   1   72  10000
'''
sum = pd.DataFrame([(sum(data["Nh"]), sum(data["Rh"]), sum(data["nh"]), sum(data["NhRh"]))], columns=["Nh","Rh","nh","NhRh"], index=["Sum"])
'''
      Nh  Rh   nh    NhRh
0  90000   4  900  125000
'''
data = data.append(sum, ignore_index=False)
'''
               Nh  Rh   nh    NhRh
Houses      35000   2  504   70000
Apartments  45000   1  324   45000
Condos      10000   1   72   10000
0           90000   4  900  125000
'''
V_p_str = (data["Nh"][0]/data["Nh"][3])**2 * 0.45*(1-.45)/350 + (data["Nh"][1]/data["Nh"][3])**2 * 0.25*(1-.25)/450 + (data["Nh"][2]/data["Nh"][3])**2 * 0.03*(1-.03)/100
# 0.000214
V_p_srs = (35/90*.45 + 45/90*.25 + 10/90*.03) * (1-(35/90*.45 + 45/90*.25 + 10/90*.03)) / 900
# 0.000234
