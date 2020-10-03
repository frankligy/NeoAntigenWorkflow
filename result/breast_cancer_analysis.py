import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc
from matplotlib.colors import ListedColormap
import collections
import math


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


## after talking with Anu on Sep 15th
# task 1: multivariable correlation in 3D plot with multiple correlation
def multiple_correlation(xy,xz,yz,n,k=2):
    # refer to http://www.real-statistics.com/correlation/multiple-correlation/
    # z is dependent variable and y,z are independent variable
    Rxyz = math.sqrt((abs(xz ** 2) + abs(yz ** 2) - 2 * xz * yz * xy) / (1 - abs(xy ** 2)))
    R2 = Rxyz ** 2
    k = 2  # Number of independent variables
    R2_adj = 1 - (((1 - R2) * (n - 1)) / (n - k - 1))
    return Rxyz,R2,R2_adj

df3d = pd.read_csv('/Users/ligk2e/Desktop/final2.txt',sep='\t',header=None)
df3d_corr = df3d.corr()
'''
1. neojunctions
2. neoantigens
3. neojunctions that give rise to neoantigens

here 2 is z, 1,3 are x and y
'''
xy = df3d_corr.loc[1,3]
xz = df3d_corr.loc[1,2]
yz = df3d_corr.loc[2,3]
metric = multiple_correlation(xy,xz,yz,len(df3d))

from mpl_toolkits import mplot3d
df3d = pd.read_csv('/Users/ligk2e/Desktop/final2.txt',sep='\t',header=None)
ax = plt.axes(projection='3d')
ax.scatter3D(df3d[1],df3d[3],df3d[2])
ax.set_xlabel('Neojunction count per patient')
ax.set_ylabel('Neojunctions that give rise to neoantigens')
ax.set_zlabel('Neoantigen count')
ax.legend(['Adjusted multiple correlation: {0:.2f}'.format(metric[2])])

# color code PAM-50 groups

# overlay basal-like, so all basal like sample will be labelled as 1
ax = plt.axes(projection='3d')
label = np.array([1 if pam50_dic[item]=='Basal-like' else 0 for item in df3d[0]])
colors = ListedColormap(['g','r'])
scatter = ax.scatter3D(df3d[1],df3d[2],df3d[3],c=label,cmap=colors)
legend1 = ax.legend(handles=scatter.legend_elements()[0],labels=['non-basal','basal'])
plt.gca().add_artist(legend1)
#ax.legend(['Multiple correlation: {0:.2f}'.format(metric[2])],loc='upper left')
plt.title('Distribution of basal patients')

# overlay HER2
ax = plt.axes(projection='3d')
label = np.array([1 if pam50_dic[item]=='HER2-enriched' else 0 for item in df3d[0]])
colors = ListedColormap(['g','r'])
scatter = ax.scatter3D(df3d[1],df3d[2],df3d[3],c=label,cmap=colors)
legend1 = ax.legend(handles=scatter.legend_elements()[0],labels=['non-HER2','HER2'])
plt.gca().add_artist(legend1)
#ax.legend(['Multiple correlation: {0:.2f}'.format(metric[2])],loc='upper left')
plt.title('Distribution of HER2 patients')

# overlay luminalA
ax = plt.axes(projection='3d')
label = np.array([1 if pam50_dic[item]=='Luminal A' else 0 for item in df3d[0]])
colors = ListedColormap(['g','r'])
scatter = ax.scatter3D(df3d[1],df3d[2],df3d[3],c=label,cmap=colors)
legend1 = ax.legend(handles=scatter.legend_elements()[0],labels=['non-LuminalA','LuminalA'])
plt.gca().add_artist(legend1)
#ax.legend(['Multiple correlation: {0:.2f}'.format(metric[2])],loc='upper left')
plt.title('Distribution of LuminalA patients')

# overlay luminal B
ax = plt.axes(projection='3d')
label = np.array([1 if pam50_dic[item]=='Luminal B' else 0 for item in df3d[0]])
colors = ListedColormap(['g','r'])
scatter = ax.scatter3D(df3d[1],df3d[2],df3d[3],c=label,cmap=colors)
legend1 = ax.legend(handles=scatter.legend_elements()[0],labels=['non-LuminalB','LuminalB'])
plt.gca().add_artist(legend1)
#ax.legend(['Multiple correlation: {0:.2f}'.format(metric[2])],loc='upper left')
plt.title('Distribution of LuminalB patients')

# overlay normal-like
ax = plt.axes(projection='3d')
label = np.array([1 if pam50_dic[item]=='Normal-like' else 0 for item in df3d[0]])
colors = ListedColormap(['g','r'])
scatter = ax.scatter3D(df3d[1],df3d[2],df3d[3],c=label,cmap=colors)
legend1 = ax.legend(handles=scatter.legend_elements()[0],labels=['non-Normal-like','Normal-like'])
plt.gca().add_artist(legend1)
#ax.legend(['Multiple correlation: {0:.2f}'.format(metric[2])],loc='upper left')
plt.title('Distribution of Normal-like patients')


## 2d with correlation
plt.scatter(df3d[1],df3d[2])
plt.xlabel('Neojunction counts')
plt.ylabel('Neoantigen counts')
corr = sc.pearsonr(df3d[1],df3d[2])
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])])
plt.title('Correlation between Neojunction count and Neoantigen count')

plt.scatter(df3d[1],df3d[3])
plt.xlabel('Neojunction counts')
plt.ylabel('Neojunctions that give rise to neoantigens')
corr = sc.pearsonr(df3d[1],df3d[3])
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])])
plt.title('Correlation between Neojunction count and Neojunctions that give rise to nenantigens')

plt.scatter(df3d[3],df3d[2])
plt.xlabel('Neojunctions that give rise to neoantigens')
plt.ylabel('Neoantigens counts')
corr = sc.pearsonr(df3d[2],df3d[3])
plt.legend(['Pearson correlation: {0:.2f}'.format(corr[0])])
plt.title('Correlation between Neojunctions that give rise to nenantigens and neoantigen counts')


#### enrichment analysis input file preparation
# matrix file
df = pd.read_csv('/Users/ligk2e/Desktop/stratification.txt',sep='\t')

df_new = df.loc[df['label']!='outliers']
df_new = df_new.set_index(pd.Index(np.arange(df_new.shape[0])))

matrix = pd.crosstab(df_new['TCGA-ID'],df_new['label'])
matrix.reset_index(inplace=True)
matrix.rename(columns={'TCGA-ID':'uid'},inplace=True)
matrix.to_csv('/Users/ligk2e/Desktop/matrix.txt',sep='\t',index=None)

# ref file
data = np.loadtxt('/Users/ligk2e/Desktop/Anu/ref1.txt').astype(np.int)
event = pd.read_csv('/Users/ligk2e/Desktop/Anu/eve.txt',sep='\t',header=None)[0]
event.name = 'event'
sample = pd.read_csv('/Users/ligk2e/Desktop/Anu/tmp2.txt',sep='\t',header=None)[0]
sample.name = 'sample'
df = pd.DataFrame(data=data,index=event,columns=sample)
df.to_csv('/Users/ligk2e/Desktop/Anu/crosstab_reference.txt',sep='\t')
df_stack = df.stack()   # a multiindex series
df_stack = df_stack.loc[df_stack==1]    # combinations that are true
new = df_stack.index.to_frame()
new['duplicate'] = new['event']
new2 = new[['sample','event','duplicate']]
new2.to_csv('/Users/ligk2e/Desktop/Anu/reference.txt',sep='\t',header=None,index=None)


## ref4000 files
data = np.loadtxt('/Users/ligk2e/Desktop/Anu/e4000/ref4000.txt').astype(np.int)
event = pd.read_csv('/Users/ligk2e/Desktop/Anu/e4000/e4000.txt',sep='\t',header=None)[0]
event.name = 'event'
sample = pd.read_csv('/Users/ligk2e/Desktop/Anu/tmp2.txt',sep='\t',header=None)[0]
sample.name = 'sample'
df = pd.DataFrame(data=data,index=event,columns=sample)
df.to_csv('/Users/ligk2e/Desktop/Anu/e4000/crosstab_reference.txt',sep='\t')
df_stack = df.stack()   # a multiindex series
df_stack = df_stack.loc[df_stack==1]    # combinations that are true
new = df_stack.index.to_frame()
new['duplicate'] = new['event']
new2 = new[['sample','event','duplicate']]
new2.to_csv('/Users/ligk2e/Desktop/Anu/e4000/reference.txt',sep='\t',header=None,index=None)


# a chi-square analysis, does pam-50 defined groups are independnt regarding their distribution in high and low group
group = pd.read_csv('/Users/ligk2e/Desktop/group.txt',sep='\t')

pam50_dic = {}
for i in range(group.shape[0]):
    subtype = group['PAM50-mRNA'].iloc[i]
    tcga = group['Complete-TCGA-ID'].iloc[i]
    pam50_dic[tcga] = subtype

df = pd.read_csv('/Users/ligk2e/Desktop/final.txt',sep='\t',header=None)

# the contigency table would be 2*4
'''
------------------------------------------------
       | basal  | her2  | luminalA | luminal B |
------------------------------------------------
high   |   20   |   7   |     30   |     32    |
------------------------------------------------
low    |   25   |  14   |   69     |      21   |
----------------------------------------------
'''

stratification = pd.read_csv('/Users/ligk2e/Desktop/stratification.txt',sep='\t')

# reverse pam50_dic
pam50_dic_reverse = {}
for k,v in pam50_dic.items():
    try:
        pam50_dic_reverse[v].append(k)
    except KeyError:
        pam50_dic_reverse[v] = []
        pam50_dic_reverse[v].append(k)

# construct a dic for high,low
hl_f = {}
hl_b = {}
for i in range(stratification.shape[0]):
    tcga = stratification['TCGA-ID'].iloc[i]
    belong = stratification['label'].iloc[i]
    try:
        hl_f[belong].append(tcga)
    except KeyError:
        hl_f[belong] = []
        hl_f[belong].append(tcga)
    hl_b[tcga] = belong


def extracting(pam50_dic_reverse,hl_b,subtype):
    # 'Basal-like' 'HER2-enriched' 'Luminal A' 'Luminal B'
    patients = pam50_dic_reverse[subtype]
    high_count,low_count = 0,0
    for i in patients:
        try:
            label = hl_b[i]
        except:
            continue
        if label == 'high':
            high_count += 1
        elif label == 'low':
            low_count += 1
    return high_count,low_count

high_basal, low_basal = extracting(pam50_dic_reverse,hl_b,'Basal-like')
high_her2, low_her2 = extracting(pam50_dic_reverse,hl_b,'HER2-enriched')
high_la,low_la = extracting(pam50_dic_reverse,hl_b,'Luminal A')
high_lb, low_lb = extracting(pam50_dic_reverse,hl_b,'Luminal B')

contigency = np.empty([2,4])
contigency[0,0] = high_basal
contigency[1,0] = low_basal
contigency[0,1] = high_her2
contigency[1,1] = low_her2
contigency[0,2] = high_la
contigency[1,2] = low_la
contigency[0,3] = high_lb
contigency[1,3] = low_lb
contigency = contigency.astype(np.int)

import scipy.stats as sc
obs_all = sc.chi2_contingency(contigency)   # for all   # 0.003 !!!!!
obs_basal_her2 = sc.chi2_contingency(contigency[:,0:2])   # 0.55
obs_basal_luminalA = sc.chi2_contingency(contigency[:,[0,2]])  # 0.14
obs_basal_luminalB = sc.chi2_contingency(contigency[:,[0,3]])  # 0.17

obs_her2_luminalA = sc.chi2_contingency(contigency[:,[1,2]])   # 0.98
obs_her2_luminalB = sc.chi2_contingency(contigency[:,[1,3]])   # 0.06 !!!!!
obs_luminalA_luminalB = sc.chi2_contingency(contigency[:,[2,3]])    # 0.00062   !!!!!!



