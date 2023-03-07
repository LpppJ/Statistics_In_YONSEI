from cmath import log
from email.mime import base
from os import stat
from secrets import choice
from statistics import correlation
from tempfile import tempdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from random import sample

# Problem 4 in Section 6.9
sample = pd.DataFrame([(7/16, 11), (3/16, 20), (3/16, 24), (3/16, 245)], columns=("psi","ti"), index=("A","B","C","D"))
sample["^t_psi"] = sample["ti"]/sample["psi"]
sample["(^t_psi -t)**2"] = (sample["^t_psi"]-np.sum(sample["ti"]))**2
sample["(^t_psi -t)**2"] = sample["(^t_psi -t)**2"].apply(lambda x: '%.4f' % x)
'''
>>> sample
      psi   ti       ^t_psi (^t_psi -t)**2
A  0.4375   11    25.142857     75546.4490
B  0.1875   20   106.666667     37377.7778
C  0.1875   24   128.000000     29584.0000
D  0.1875  245  1306.666667   1013377.7778
'''


E_t_psi = sum(sample["psi"][i]*sample["^t_psi"][i] for i in range(0,len(sample)))
'''
>>> E_t_psi
300.0
'''

V_t_psi = sum(sample["psi"][i]*float(sample["(^t_psi -t)**2"][i]) for i in range(0,len(sample)))
'''
>>> V_t_psi
235615.2381125
'''
# Conclusion :
# Psi = (7/16,3/16,3/16,3/16) produces a higher variance than simple random sampling


# Problem 13 in Section 6.9
# (a)
statepop = pd.read_csv("statepop.csv")
statepop["psi"] = statepop["popn"]/sum(statepop["popn"])
plt.plot(statepop["veterans"], statepop["psi"])
plt.show()
np.corrcoef(statepop["veterans"], statepop["popn"]) # 0.98719407
'''
"veterns vs psi" looks like a 45-degree line.
and corr(veterans, popn) = 0.98719407 is closed to 1
Thus, we expect unequal probability sampling to be efficient.
'''

# (b)
statepop["veterans/psi"] = statepop["veterans"]/statepop["psi"]

t_veterans_hat_psi = sum(statepop["veterans/psi"])/len(statepop) # 11476938.645127924

SE_t_hat_psi = sum(np.sqrt((statepop["veterans/psi"][i]-t_veterans_hat_psi)**2/(len(statepop)-1)/len(statepop))
                   for i in range(0,len(statepop))) # 2071043.4012261478

# (c)
statepop["Viet_veterans"] = statepop["veterans"]*statepop["percviet"]/100
statepop["Viet_veterans/psi"] = statepop["Viet_veterans"]/statepop["psi"]

t_Viet_veterans_hat_psi = sum(statepop["Viet_veterans/psi"])/len(statepop) #3309960.37216951

SE_t_Viet_veterans_hat_psi = sum(np.sqrt((statepop["Viet_veterans/psi"][i]-t_Viet_veterans_hat_psi)**2/(len(statepop)-1)/len(statepop)) for i in range(0,len(statepop))) # 694614.0313025394

# Problem 36 in Section 5.8
baseball = pd.read_csv("baseball.csv",names=("team","leagueID","player","salary","POS","G","GS","InnOuts","PO","A","E","DP","PB","GF","AB","R","H","SecB","ThiB","HR","RBI","SB","CS","BB","SO","IBB","HBP","SH","SF","GIDP"))
baseball["logsal"] = np.log(baseball["salary"])



# (a)
'''
In the baseball.csv file, there are data of 797 players from 30 teams. Since the team is psu and contains data of less than 30 players in a team, you have to choose at least 6 teams to use the data of 150 players. 6 teams will be selected by SRS.
'''



# (b)
team_list =[]
for name in baseball["team"]:
    if name not in team_list:
        team_list.append(name)
team_list
'''
['ANA', 'ARI', 'ATL', 'BAL', 'BOS', 'CHA', 'CHN', 'CIN', 'CLE', 'COL', 'DET', 'FLO', 'HOU', 'KCA', 'LAN', 'MIL', 'MIN', 'MON', 'NYA', 'NYN', 'OAK', 'PHI', 'PIT', 'SDN', 'SEA', 'SFN', 'SLN', 'TBA', 'TEX', 'TOR']
'''
sample_team = sample(team_list, 6) # Random result : ['FLO', 'CHA', 'TBA', 'SDN', 'PIT', 'ANA']
sample_data = baseball[baseball["team"].isin(sample_team)]

'''
>>> sample_data
    team leagueID    player    salary POS    G   GS  InnOuts   PO  ...  CS  BB  SO  IBB   HBP  SH   SF  GIDP     logsal
0    ANA       AL  anderga0   6200000  CF  112   92     2375  211  ...   1  29  75  6.0   1.0   0  3.0     3  15.640060
1    ANA       AL  colonba0  11000000   P    3   34      625    8  ...   0   0   1  0.0   0.0   0  0.0     0  16.213406 
..   ...      ...       ...       ...  ..  ...  ...      ...  ...  ...  ..  ..  ..  ...   ...  ..  ...   ...        ...  
742  TBA       AL  standja0    302500   P    3    1       30    2  ...   0   0   0  0.0   0.0   0  0.0     0  12.619837  
743  TBA       AL  zambrvi0    325000   P    3   22      384   10  ...   0   0   1  0.0   0.0   0  0.0     0  12.691580  
[157 rows x 31 columns]
'''

sample_data_0 = sample_data[sample_data["team"]==sample_team[0]]
sample_data_1 = sample_data[sample_data["team"]==sample_team[1]]
sample_data_2 = sample_data[sample_data["team"]==sample_team[2]]
sample_data_3 = sample_data[sample_data["team"]==sample_team[3]]
sample_data_4 = sample_data[sample_data["team"]==sample_team[4]]
sample_data_5 = sample_data[sample_data["team"]==sample_team[5]]

plt.boxplot([sample_data_0["logsal"], sample_data_1["logsal"], sample_data_2["logsal"], sample_data_3["logsal"], sample_data_4["logsal"], sample_data_5["logsal"],])
plt.xticks([1,2,3,4,5,6],[sample_team[0], sample_team[1], sample_team[2], sample_team[3], sample_team[4], sample_team[5]])
plt.show()



# (c)
# We will use 'Unbiased estimation'
N=30
n=6
M0 = len(baseball) # 797
t_hat_unb = np.sum(sample_data["logsal"]) * N/n # 10842.9364
y_bar_hat_unb = t_hat_unb / M0 # 13.6046

t_hat=[]
t_hat.append(np.sum(sample_data_0["logsal"])) # 356.9215 : The sum of the logsal of the players whose team is FlO(=sample_team[0])
t_hat.append(np.sum(sample_data_1["logsal"])) # 362.2564 : The sum of the logsal of the players whose team is CHA(=sample_team[1])
t_hat.append(np.sum(sample_data_2["logsal"])) # 350.9497 : The sum of the logsal of the players whose team is TBA(=sample_team[2])
t_hat.append(np.sum(sample_data_3["logsal"])) # 362.6829 : The sum of the logsal of the players whose team is SDN(=sample_team[3])
t_hat.append(np.sum(sample_data_4["logsal"])) # 361.3865 : The sum of the logsal of the players whose team is PIT(=sample_team[4])
t_hat.append(np.sum(sample_data_5["logsal"])) # 374.3900 : The sum of the logsal of the players whose team is ANA(=sample_team[5])
# t_hat = [356.9215, 362.2564, 350.9497, 362.6829, 361.3865, 374.3900]

s_t_square = sum((t_hat[i]-t_hat_unb/N)**2 for i in range(0,len(sample_team)))/(n-1) # 60.0760
SEhat_t_hat_unb = N * np.sqrt((1-n/N)*s_t_square/n) # 84.9065
SEhat_y_bar_hat_unb = SEhat_t_hat_unb/M0 # 0.1065
CI_y_bar_hat_unb = [y_bar_hat_unb-1.96*SEhat_y_bar_hat_unb, y_bar_hat_unb+1.96*SEhat_y_bar_hat_unb]
# [13.395884152239386, 13.813492234191948]



# (d)
# We will use 'Ratio estimation'
# We expect 'Number of pitchers' to be positively correlated with 'Mi'(ith psu size)
sample_pitchers = sample_data[sample_data["POS"]=="P"]
yr_bar_hat = len(sample_pitchers)/len(sample_data) # 73/157 = 0.4649
Mi = [len(sample_data_0), len(sample_data_1), len(sample_data_2), len(sample_data_3), len(sample_data_4), len(sample_data_5)]
# [26, 26, 26, 26, 27, 26]

y_bar=[]
y_bar.append(len(sample_data_0[sample_data_0['POS']=='P'])/len(sample_data_0))
# 0.4615 : The ratio of the Pitcher of the players whose team is FlO(=sample_team[0])
y_bar.append(len(sample_data_1[sample_data_1['POS']=='P'])/len(sample_data_1))
# 0.5 : The ratio of the Pitcher of the players whose team is CHA(=sample_team[1])
y_bar.append(len(sample_data_2[sample_data_2['POS']=='P'])/len(sample_data_2))
# 0.4615 : The ratio of the Pitcher of the players whose team is TBA(=sample_team[2])
y_bar.append(len(sample_data_3[sample_data_3['POS']=='P'])/len(sample_data_3))
# 0.4615 : The ratio of the Pitcher of the players whose team is SDN(=sample_team[3])
y_bar.append(len(sample_data_4[sample_data_4['POS']=='P'])/len(sample_data_4))
# 0.4444 : The ratio of the Pitcher of the players whose team is PIT(=sample_team[4])
y_bar.append(len(sample_data_5[sample_data_5['POS']=='P'])/len(sample_data_5))
# 0.4615 : The ratio of the Pitcher of the players whose team is ANA(=sample_team[5])

# y_bar = [0.4615, 0.5, 0.4615, 0.4615, 0.4444, 0.4615]
SE_yr_bar_hat = sum(np.sqrt((1-n/N)/(n*np.mean(Mi)**2) * (Mi[i]**2)*((y_bar[i]-yr_bar_hat)**2)/(n-1)) for i in range(0,len(sample_team)))
# 0.011368478501131703

CI_yr_bar_hat = [yr_bar_hat-1.96*SE_yr_bar_hat, yr_bar_hat+1.96*SE_yr_bar_hat]
#[0.4426859350040239, 0.48725037072846017]



# (e)
# Excercise 32 in Chapter 2

# 2-32-(a)
baseball = pd.read_csv("baseball.csv", names=("team","leagueID","player","salary","POS","G","GS","InnOuts","PO","A","E","DP","PB","GF","AB","R","H","SecB","ThiB","HR","RBI","SB","CS","BB","SO","IBB","HBP","SH","SF","GIDP"))
sample_data_SRS = baseball.sample(n=150)

# 2-32-(b)
sample_data_SRS["logsal"] = np.log(sample_data_SRS["salary"])
plt.hist(sample_data_SRS["logsal"])
plt.show() # Not normal distribution

# 2-32-(c)
np.mean(sample_data_SRS["logsal"]) # 13.8840

s_square = sum((list(sample_data_SRS["logsal"])[i]-np.mean(sample_data_SRS["logsal"]))**2 for i in range(0,len(sample_data_SRS)))/(150-1)
# 1.4548
SE_ybar = np.sqrt((1-150/797)*s_square/150) # 0.08873

CI = [np.mean(sample_data_SRS["logsal"])-1.96*SE_ybar, np.mean(sample_data_SRS["logsal"])+1.96*SE_ybar]
# [13.710082806267966, 14.057918568244602]

# 2-32-(d)
sample_pitchers_SRS = sample_data_SRS[sample_data_SRS["POS"]=="P"]
p_hat = len(sample_pitchers_SRS)/len(sample_data_SRS) # 0.52

SE_p_hat = np.sqrt((1-150/797)*p_hat*(1-p_hat)/(150-1)) # 0.03687670727684561

CI = [p_hat-1.96*SE_p_hat, p_hat+1.96*SE_p_hat] # [0.4477216537373826, 0.5922783462626174]
'''
When comparing the estimated values derived from #36 of chapter 5 and #32 of chapter 2,
the confidence interval of cluster sampling is narrower under the same significance level.
It is because of the intergroup homogeneity of the sample.
'''
