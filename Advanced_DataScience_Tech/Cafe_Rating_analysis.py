import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyparsing import line
from sklearn import linear_model
import statsmodels.api as sm
from random import randint
from PIL import Image

# Problem 1
cafe = pd.read_csv("cafe.csv", encoding="cp949")
user_review = pd.read_csv("user_reviews.csv", encoding="cp949")
user_review.groupby('place_id').size()
'''
place_id
11853424       93
13455155      102
13546171      204
20108351      153
20620536      212
             ...
1958401576    549
1969214801    462
1980225181    103
1987513835    132
1989383767     96
Length: 124, dtype: int64
'''
index = np.random.choice(np.arange(len(cafe)))
mychoice = [cafe['place_id'][index], cafe['name'][index]]
lst = list(np.where(user_review['place_id'] == mychoice[0])[0])
total = 0
for i in lst:
  total += user_review['user_avg_rating'][i]

mychoice[1]
ref = round(total/len(lst),4)

# In this attempt,  
mychoice[1] # '매머드익스프레스시청점'
ref # 4.384

myImage = Image.open("HW3_Q1.jpg")
myImage.show()
# The above figure is the process of deriving the following formula.

lambda_ = 0.01
x_a = cafe['avg_rating']
n_a = len(x_a)
for i in range(100):
    g = 1/n_a * sum(np.exp(-1/n_a * lambda_ * x_a) * x_a) - ref
    der_g = 1/n_a * sum(np.exp(-1/n_a * lambda_ * x_a) * (-1/n_a * x_a) * x_a)
    
    step = - g / der_g
    lambda_ = lambda_ + step
lambda_ # 0.5312
step # 0.0

w = np.exp(-1/n_a * lambda_ * x_a)
adj_w = w/sum(w) * n_a
adj_w
cafe['adj_rating'] = adj_w*x_a
cafe.sort_values('adj_rating', ascending=False)
'''
>>> cafe.sort_values('adj_rating', ascending=False)
       place_id        name        ...  user_avg_rating   adj_rating
55   1591980712    커피를그리다    ...     4.904469         4.942262       
92   1169675907  홀빈커피작은가게  ...     4.667857         4.913421     
73     33877615       아필립       ...     4.698374         4.889368
57   1608809718      에그스팟      ...     4.849515         4.884315
101  1325051960      카페스타      ...     4.783367         4.871728
..          ...         ...        ...       ...               ...
50   1267517766      힐링커피      ...     4.148000         4.175342
8      35516328   아재커피시청점   ...     4.283893         4.161613      
0      33434383   커피마마종암점   ...     4.053279         4.161316      
100  1324604066  마이브커피이야기  ...     4.077143         4.067820     
102  1332877877      이트커피  ...        3.959701          3.974965
[124 rows x 9 columns]
'''

# Problem 2
# 1. Gradient Descent Algorithm

# (STEP1) Checking for Linearity
Data = pd.read_csv("3rd_example.csv")
y=Data["y"]
z1=Data["z1"]
z2=Data["z2"]
z=Data.loc[:,["z1","z2"]]
plt.scatter(y,z1,color="orange",label="y vs z1")
plt.scatter(y,z2,color="green",label="y vs z2")
plt.legend
plt.show()

# (STEP2) Initial values : Beta, Learning rate
Beta = [0, 0, 0] # [beta0, beta1, beta2]

L = 0.001 # The Learning rate
n = int(len(z1)) # 400

# (STEP3) Initial Iteration
Y_pred = Beta[0] + Beta[1]*z1 + Beta[2]*z2
D_Beta0 = (-2/n) * sum(y-Y_pred)
D_Beta1 = (-2/n) * sum(z1*(y-Y_pred))
D_Beta2 = (-2/n) * sum(z2*(y-Y_pred))
Beta_new = [Beta[0]-L*D_Beta0,  Beta[1]-L*D_Beta1, Beta[2]-L*D_Beta2]
# [0.003317, 0.006221, 0.007350]
Y_pred_new = Beta_new[0] + Beta_new[1]*z1 + Beta_new[2]*z2


# (STEP4) Second ~ Last Iteration
while np.sqrt(sum((Y_pred_new[i]-Y_pred[i])**2 for i in range(0,int(n)))) >= 1e-06:
# while "Distance between Y_pred_new and Y_pred" is greater than or equal to 1e-06

    # Updating The Learing rate
    alpha = np.random.uniform(0,0.5)
    beta = np.random.uniform(0,1)
    while 1/n*sum((y-Y_pred_new)**2) > 1/n*sum((y-Y_pred)**2) - alpha*L*(D_Beta0**2+D_Beta1**2+D_Beta2**2):
        L = L*beta
    
    Beta = Beta_new    
    Y_pred = Y_pred_new
    D_Beta0 = (-2/n) * sum(y-Y_pred)
    D_Beta1 = (-2/n) * sum(z1*(y-Y_pred))
    D_Beta2 = (-2/n) * sum(z2*(y-Y_pred))
    Beta_new[0] = Beta[0] - L* D_Beta0
    Beta_new[1]= Beta[1] - L* D_Beta1
    Beta_new[2] = Beta[2] - L* D_Beta2
    Y_pred_new = Beta_new[0] + Beta_new[1]*z1 + Beta_new[2]*z2
    
print(Beta_new)
# beta0 = 0.5877
# beta1 = 0.4197
# beta2 = 0.4462

# (STEP5) Checking the estimator of Beta
model = sm.OLS(y, sm.add_constant(z))
result = model.fit()
result.summary()
'''
                            OLS Regression Results
==============================================================================
Dep. Variable:                      y   R-squared:                       0.423
Model:                            OLS   Adj. R-squared:                  0.420
Method:                 Least Squares   F-statistic:                     145.4
Date:                Wed, 15 Jun 2022   Prob (F-statistic):           4.28e-48
Time:                        12:07:12   Log-Likelihood:                -661.21
No. Observations:                 400   AIC:                             1328.
Df Residuals:                     397   BIC:                             1340.
Df Model:                           2
Covariance Type:            nonrobust
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
const          0.5878      0.095      6.207      0.000       0.402       0.774
z1             0.4197      0.046      9.195      0.000       0.330       0.509
z2             0.4462      0.054      8.319      0.000       0.341       0.552
==============================================================================
Omnibus:                        3.276   Durbin-Watson:                   2.102
Prob(Omnibus):                  0.194   Jarque-Bera (JB):                3.231
Skew:                           0.179   Prob(JB):                        0.199
Kurtosis:                       2.744   Cond. No.                         4.14
==============================================================================
'''
# (STEP6) Conclusion : In OLS Regression Results, we can find "coef of [const/z1/z2]", which are estimated by "[beta0/beta1/beta2]"

# 2. Stochastic Gradient Descent Algorithm

# (STEP1) Initial values : Beta, Learning rate
Beta = [0, 0, 0] # [beta0, beta1, beta2]

L = 0.001 # The Learning rate
n = int(len(z)) # 400

# (STEP2) Initial Iteration
Y_pred = Beta[0] + Beta[1]*z1 + Beta[2]*z2
D_Beta0 = (-2/n) * sum(y-Y_pred)
D_Beta1 = (-2/n) * sum(z1*(y-Y_pred))
D_Beta2 = (-2/n) * sum(z2*(y-Y_pred))
Beta_new = [Beta[0]-L*D_Beta0,  Beta[1]-L*D_Beta1, Beta[2]-L*D_Beta2]
# [0.003317, 0.006221, 0.007350]
Y_pred_new = Beta_new[0] + Beta_new[1]*z1 + Beta_new[2]*z2

# (STEP3) Second ~ Last Iteration
while np.sqrt(sum((Y_pred_new[k]-Y_pred[k])**2 for k in range(0,int(n)))) >= 1e-06:
# while "Distance between Y_pred_new and Y_pred" is greater than or equal to 1e-06

    # Updating The Learing rate
    alpha = np.random.uniform(0,0.5)
    beta = np.random.uniform(0,1)
    while 1/n*sum((y-Y_pred_new)**2) > 1/n*sum((y-Y_pred)**2) - alpha*L*(D_Beta0**2+D_Beta1**2+D_Beta2**2):
        L = L*beta
    
    Beta = Beta_new    
    Y_pred = Y_pred_new
    i = np.random.randint(0,n)
    D_Beta0 = (-2/n) * n* (y[i]-Y_pred[i])
    D_Beta1 = (-2/n) * n* z1[i]*(y[i]-Y_pred[i])
    D_Beta2 = (-2/n) * n* z2[i]*(y[i]-Y_pred[i])
    Beta_new[0] = Beta[0] - L*D_Beta0
    Beta_new[1]= Beta[1] - L*D_Beta1
    Beta_new[2] = Beta[2] - L*D_Beta2
    Y_pred_new = Beta_new[0] + Beta_new[1]*z1 + Beta_new[2]*z2
    
print(Beta_new)
# beta0 = 0.58727
# beta1 = 0.41969
# beta2 = 0.44643

# (STEP4) Conclusion : We can estimate "[beta0/beta1/beta2]" which are also near by "coef of [const/z1/z2]"

# Problem 3
model = sm.OLS(z2, sm.add_constant(z1))
result = model.fit()
result.summary()
"""
                            OLS Regression Results

==============================================================================
Dep. Variable:                     z2   R-squared:
  0.223
Model:                            OLS   Adj. R-squared:
  0.221
Method:                 Least Squares   F-statistic:
  114.1
Date:                Fri, 17 Jun 2022   Prob (F-statistic):           1.38e-23
Time:                        23:45:06   Log-Likelihood:
-634.60
No. Observations:                 400   AIC:
  1273.
Df Residuals:                     398   BIC:
  1281.
Df Model:                           1

Covariance Type:            nonrobust

==============================================================================
                 coef    std err          t      P>|t|      [0.025     
 0.975]
------------------------------------------------------------------------------
const          1.0865      0.070     15.575      0.000       0.949     
  1.224
z1             0.4017      0.038     10.682      0.000       0.328       0.476
==============================================================================
Omnibus:                        3.725   Durbin-Watson:                   2.099
Prob(Omnibus):                  0.155   Jarque-Bera (JB):                3.600
Skew:                           0.186   Prob(JB):                        0.165
Kurtosis:                       2.720   Cond. No.                         2.40
==============================================================================

Notes:
[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
"""

model = sm.OLS(y, sm.add_constant(z1))
result = model.fit()
result.summary()
"""
                            OLS Regression Results
==============================================================================
Dep. Variable:                      y   R-squared:                       0.322
Model:                            OLS   Adj. R-squared:                  0.320
Method:                 Least Squares   F-statistic:                     189.1
Date:                Fri, 17 Jun 2022   Prob (F-statistic):           1.75e-35
Time:                        23:45:54   Log-Likelihood:                -693.34
No. Observations:                 400   AIC:                             1391.
Df Residuals:                     398   BIC:                             1399.
Df Model:                           1
Covariance Type:            nonrobust
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
const          1.0726      0.081     13.276      0.000       0.914       1.231
z1             0.5989      0.044     13.752      0.000       0.513       0.685
==============================================================================
Omnibus:                        1.354   Durbin-Watson:                   2.120
Prob(Omnibus):                  0.508   Jarque-Bera (JB):                1.317
Skew:                           0.140   Prob(JB):                        0.518
Kurtosis:                       2.980   Cond. No.                         2.40
==============================================================================

Notes:
[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
"""
# The above figure is the estimation process of theta
