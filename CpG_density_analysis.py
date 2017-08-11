
# coding: utf-8

# In[6]:

import matplotlib.pyplot as plt
import numpy as np
get_ipython().magic('matplotlib inline')


# In[1]:

import pickle 
cancerGeneToCpGDensity = pickle.load(open('/Users/khandekara2/Documents/methylationProject/01_data/cancerGeneToCpGDensity.pickle', 'rb'))
noncancerGeneToCpGDensity = pickle.load(open('/Users/khandekara2/Documents/methylationProject/01_data/noncancerGeneToCpGDensity.pickle', 'rb'))


# In[2]:

print (cancerGeneToCpGDensity)


# In[4]:

cancer_ratios = [i * 100 for i in cancerGeneToCpGDensity.values()]
noncancer_ratios = [i * 100 for i in noncancerGeneToCpGDensity.values()]


# In[10]:

plt.hist(cancer_ratios, bins=12)


# In[9]:

plt.hist(noncancer_ratios, bins=25)


# In[ ]:



