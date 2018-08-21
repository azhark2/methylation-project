import numpy as np
import pandas as pd
import pickle
from scipy.stats import chi2_contingency
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import r2_score
from scipy.stats.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import sys
import decimal
from scipy import stats
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# slope, intercept, r_value, p_value, std_err = stats.linregress(vafs, drop)
# print (slope)
# print (intercept)
# print (r_value ** 2)
# print (p_value)
# print (std_err)

#create first plot of vaf vs. drop
mutated = pd.read_csv('/data/khandekara2/imputation/mutations_reconstructed.tsv', sep='\t')
m = mutated[mutated['drop'] >= 0]
g = (sns.jointplot('vaf', # Horizontal axis
           'drop', # Vertical axis
           data=m, # Data source
           size=4)).set_axis_labels("Variant Allele Frequency", "Drop in Methylation Ratio")
g.ax_joint.cla()
plt.sca(g.ax_joint)
plt.scatter(m.vaf, m['drop'], s=10)
plt.title('')
plt.xlabel('Variant Allele Frequency')
plt.ylabel('Drop in Methylation Ratio')
g.savefig('vaf_vs_drop.png', dpi=600, bbox_inches='tight')

#create second plot of prior ratio vs. actual ratio
g = (sns.jointplot('methylation_ratio', # Horizontal axis
           'prior_ratio', # Vertical axis
           data=mutated, # Data source
           size=4)).set_axis_labels("Prior Methylation Ratio\n(Predicted)", "Observed Ratio after Mutation") # S marker size

g.ax_joint.cla()
plt.sca(g.ax_joint)
plt.scatter(m.methylation_ratio, m.prior_ratio, s=10)

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax = plt.gca()
# now plot both limits against each other
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
plt.title('')
plt.xlabel("Prior Methylation Ratio\n(Predicted)")
plt.ylabel("Observed Ratio after Mutation")
plt.savefig('actual_vs_prior.png', dpi=600)
