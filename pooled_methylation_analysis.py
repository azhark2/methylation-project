
# coding: utf-8

# In[31]:

import pickle 
WGBS_avgRatio = pickle.load(open('/Users/khandekara2/Documents/methylationProject/01_data/WGBS_avgRatio.pickle', 'rb'), encoding='latin1')
all_SNPS = pickle.load(open('/Users/khandekara2/Documents/methylationProject/01_data/all_SNPS.pickle', 'rb'), encoding='latin1')


# In[32]:

std_devs = []
for value in WGBS_avgRatio.values():
    if not np.isnan(value[1]):
        std_devs.append(value[1])
print (len(WGBS_avgRatio.keys()))


# In[34]:

import matplotlib.pyplot as plt
import numpy as np
get_ipython().magic('matplotlib inline')
plt.hist(std_devs)
plt.title('Distribution of Standard Deviations of Methylation Ratios')
plt.xlabel('Standard Deviation')
plt.ylabel('Frequency')


# In[35]:

std_devs_no_SNPS = []
for key, value in WGBS_avgRatio.items():
    key = (key[0], str(key[1]), str(key[2]))
    if key not in all_SNPS:
        if not np.isnan(value[1]):
            std_devs_no_SNPS.append(value[1])

print (len(WGBS_avgRatio.keys()))
print (len(std_devs_no_SNPS))
print (len(std_devs))


# In[36]:

print (type(key))
print (type(key[0]))
print (type(key[1]))
print (type(key[2]))


# In[37]:

loc = all_SNPS.pop()


# In[38]:

print (type(loc))
print (type(loc[0]))
print (type(loc[1]))
print (type(loc[2]))


# In[39]:

print (key)
print (key[0])
print (key[1])
print (key[2])


# In[40]:

print (loc)
print (loc[0])
print (loc[1])
print (loc[2])


# In[41]:

import matplotlib.pyplot as plt
import numpy as np
get_ipython().magic('matplotlib inline')
plt.hist(std_devs_no_SNPS)
plt.title('Distribution of Standard Deviations of Methylation Ratios with SNPS and Mutations Removed')
plt.xlabel('Standard Deviation')
plt.ylabel('Frequency')


# In[ ]:



