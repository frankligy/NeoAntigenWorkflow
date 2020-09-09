import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc
from matplotlib.colors import ListedColormap
import collections

df = pd.read_csv('/Users/ligk2e/Desktop/final.txt',sep='\t',header=None)
group = pd.read_csv('/Users/ligk2e/Desktop/group.txt',sep='\t')

plt.hist(df[1])
plt.title('Neojunction count across 804 TCGA breast cancer samples')
plt.xlabel('Neojunction counts')

plt.hist(df[2])
plt.title('AS-derived Neoantigen count across 804 TCGA breast cancer samples')
plt.xlabel('Neoantigen counts')



# correlation
plt.scatter(df[1],df[2])
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df[1],df[2])
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])])
plt.title('Correlation between Neojunction count and Neoantigen count')

# let's inspect where PAM50 clusters locate
pam50_dic = {}
for i in range(group.shape[0]):
    subtype = group['PAM50-mRNA'].iloc[i]
    tcga = group['Complete-TCGA-ID'].iloc[i]
    pam50_dic[tcga] = subtype

# overlay basal-like, so all basal like sample will be labelled as 1
label = np.array([1 if pam50_dic[item]=='Basal-like' else 0 for item in df[0]])
colors = ListedColormap(['g','r'])
scatter = plt.scatter(df[1],df[2],c=label,cmap=colors)
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df[1],df[2])
legend1 = plt.legend(handles=scatter.legend_elements()[0],labels=('non-basal-like','basal-like'))
plt.gca().add_artist(legend1)
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])],loc='lower right')
plt.title('Distribution of basal-like patients')

# overlay her-2
label = np.array([1 if pam50_dic[item]=='HER2-enriched' else 0 for item in df[0]])
colors = ListedColormap(['g','r'])
scatter = plt.scatter(df[1],df[2],c=label,cmap=colors)
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df[1],df[2])
legend1 = plt.legend(handles=scatter.legend_elements()[0],labels=('non-HER2','HER2'))
plt.gca().add_artist(legend1)
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])],loc='lower right')
plt.title('Distribution of HER2 patients')

# overlay luminalA

label = np.array([1 if pam50_dic[item]=='Luminal A' else 0 for item in df[0]])
colors = ListedColormap(['g','r'])
scatter = plt.scatter(df[1],df[2],c=label,cmap=colors)
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df[1],df[2])
legend1 = plt.legend(handles=scatter.legend_elements()[0],labels=('non-LuminalA','LuminalA'))
plt.gca().add_artist(legend1)
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])],loc='lower right')
plt.title('Distribution of LuminalA patients')

# overlay luminalB
label = np.array([1 if pam50_dic[item]=='Luminal B' else 0 for item in df[0]])
colors = ListedColormap(['g','r'])
scatter = plt.scatter(df[1],df[2],c=label,cmap=colors)
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df[1],df[2])
legend1 = plt.legend(handles=scatter.legend_elements()[0],labels=('non-LuminalB','LuminalB'))
plt.gca().add_artist(legend1)
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])],loc='lower right')
plt.title('Distribution of LuminalB patients')

# overlay normal-like
label = np.array([1 if pam50_dic[item]=='Normal-like' else 0 for item in df[0]])
colors = ListedColormap(['g','r'])
scatter = plt.scatter(df[1],df[2],c=label,cmap=colors)
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df[1],df[2])
legend1 = plt.legend(handles=scatter.legend_elements()[0],labels=('non-normal','normal-like'))
plt.gca().add_artist(legend1)
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])],loc='lower right')
plt.title('Distribution of Normal-like patients')


# most recurrent events
events = pd.read_csv('/Users/ligk2e/Desktop/all_events.txt')
cond = [False if item == 'UID' else True for item in events['UID']]
count = collections.Counter(events['UID'].loc[cond])
n50 = count.most_common(50)
value = [int(item[1]) for item in n50]
key = [item[0] for item in n50]
key_foreground = [item.split('|')[0] for item in key]
key_foreground = [item.split(':')[0]+':'+':'.join(item.split(':')[2:])for item in key_foreground]

x = np.arange(50,0,-1)
plt.barh(x,value)
plt.yticks(np.arange(50,0,-1),key_foreground,fontsize=6)
plt.title('Most recurrent tumor-specific splicing event')
plt.xlabel('Occurence')


# stratify patients
df_new = df.sort_values(2,ascending=False)
df_new = df_new.set_index(pd.Index(np.arange(df_new.shape[0])))
high = df_new.iloc[0:268]
medium = df_new.iloc[268:536]
low = df_new.iloc[536:]


