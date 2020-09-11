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
high.to_csv('/Users/ligk2e/Desktop/high.txt',sep='\t',index=None)
medium.to_csv('/Users/ligk2e/Desktop/medium.txt',sep='\t',index=None)
low.to_csv('/Users/ligk2e/Desktop/low.txt',sep='\t',index=None)

# AFTER talking with Anu on Sep 10th

# draw 3D scatter plot
from mpl_toolkits import mplot3d
df3d = pd.read_csv('/Users/ligk2e/Desktop/final2.txt',sep='\t',header=None)
ax = plt.axes(projection='3d')
ax.scatter3D(df3d[1],df3d[3],df3d[2])
ax.set_xlabel('Neojunction count per patient')
ax.set_ylabel('Neojunctions that give rise to neoantigens')
ax.set_zlabel('Neoantigen count')


# recurrent event make them unique
events = pd.read_csv('/Users/ligk2e/Desktop/all_events.txt')
events_unique = list(set(list(events['UID'].values)))   # 10270 events
cond = [False if item == 'UID' else True for item in events['UID']]
count = collections.Counter(events['UID'].loc[cond])
n500 = count.most_common(4000)
n500_array = np.array(n500)
df = pd.DataFrame({'event':n500_array[:,0],'occurence':n500_array[:,1]})
df.to_csv('/Users/ligk2e/Desktop/events4000.txt',sep='\t',index=None)


## after running on linux cluster, we get average neoantigen generations for each recurrent event
ave = pd.read_csv('/Users/ligk2e/Desktop/tmp13.txt',sep='\t',header=None)
ave.columns = ['event','occurence','average-neoantigen-generation']
ave.to_csv('/Users/ligk2e/Desktop/recurrent.txt',sep='\t',index=None)

# interquartile range and outliers
plt.boxplot(df[2])
x = np.random.normal(1,0.02,len(df[2]))
plt.plot(x,df[2],'r.',alpha=0.2)
plt.title('Neoantigen counts distribution')

q1 = df[2].quantile(q=0.25)
q3 = df[2].quantile(q=0.75)
iqr = q3-q1
upper_whisker = q3 + 1.5*iqr
lower_whisker = q1 - 1.5*iqr
# [lower_whisker, q1, median, q3, upper_whisker]
label = []
for i in range(df.shape[0]):
    if df[2].iloc[i] > lower_whisker and df[2].iloc[i] <= q1:
        label.append('low')
    elif df[2].iloc[i] > q1 and df[2].iloc[i] <= q3:
        label.append('mediem')
    elif df[2].iloc[i] > q3 and df[2].iloc[i] <= upper_whisker:
        label.append('high')
    else:
        label.append('outliers')
df['label'] = label
df.columns = ['TCGA-ID','neojunction-count','neoantigen-count','label']
df.to_csv('/Users/ligk2e/Desktop/stratification.txt',sep='\t',index=None)



