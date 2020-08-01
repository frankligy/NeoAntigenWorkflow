import numpy as np
import os
import pandas as pd
import argparse


def remove_unuseful_count(df_count):
    cond = []
    for row in df_count.itertuples():
        uid = row.UID   # row will be a namedtuple, (index=,UID=,sample=)
        if ':' in uid: cond.append(True)
        else: cond.append(False)
    df_count = df_count.loc[pd.Series(cond)]
    df_count = df_count.set_index(pd.Index(np.arange(df_count.shape[0])))
    return df_count



def construct_count_search_space(df_count):
    df_count = remove_unuseful_count(df_count)
    dic = {}
    for row in df_count.itertuples():
        ensg = row.UID.split(':')[0]   
        event = row.UID.split(':')[1:]
        event = ''.join(event)  # E2-E4 or E2-ENSG:E4
        counts = int(row[2])
        try:
            dic[ensg].append((event,counts))
        except KeyError:
            dic[ensg] = []
            dic[ensg].append((event,counts))

    return dic


def filter_zero_PSI(df_PSI):

    # change all nan to 0.0
    tumorPSI = df_PSI.iloc[:,1]
    tumorPSI[np.isnan(tumorPSI)] = 0.0
    df_PSI.iloc[:,1] = tumorPSI
  

    cond = []
    for row in df_PSI.itertuples():
        psi = row[2]
        if psi == 0: cond.append(False)
        else: cond.append(True)
    df_PSI = df_PSI[pd.Series(cond)]
    df_PSI = df_PSI.set_index(pd.Index(np.arange(df_PSI.shape[0])))

    return df_PSI
    






def main(args):
    global count
    global PSI
    global samplePSI
    global sampleCount
    global outdir

    count = args.count
    PSI = args.PSI
    samplePSI = args.samplePSI
    sampleCount = args.sampleCount
    outdir = args.outdir

    df_count = pd.read_csv(count,sep='\t',usecols=['UID',sampleCount])
    df_PSI = pd.read_csv(PSI,sep='\t',usecols=['UID',samplePSI])

    dic = construct_count_search_space(df_count)
    df_PSI = filter_zero_PSI(df_PSI)


    col = []
    for row in df_PSI.itertuples():  
        ensg = row.UID.split('|')[0].split(':')[1]
        event = row.UID.split('|')[0].split(':')[2:]
        event = ''.join(event)
        search = dic[ensg]   # narrow down search space to a certain ensg
        for item in search:
            if event == item[0]:
                read_count = item[1]
                col.append(read_count)
                break
            else:
                continue
    
    df_PSI = df_PSI.join(pd.Series(col,name='count'))
    df_PSI.to_csv(os.path.join(outdir,'{0}_ori.txt'.format(samplePSI)),sep='\t',index=None)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get UID, read count and PSI(not zero)')
    parser.add_argument('--count',type=str,default='.',help='path to your read count file, should under /AltExpression/pre-filtered/counts/')
    parser.add_argument('--PSI',type=str,default='.',help='path to your PSI EventAnnotation file, should under /AltResult/AlternativeOutput/')
    parser.add_argument('--samplePSI',type=str,default=None,help='specify which sample you want to inspect in PSI file')
    parser.add_argument('--sampleCount',type=str,default=None,help='specify which sample you want to inspect in count file')
    parser.add_argument('--outdir',type=str,default='.',help='path to the output file')
    args = parser.parse_args()
    main(args)
