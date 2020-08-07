import os
import pandas as pd
import numpy as np




if __name__ == '__main__':
    id_ = 'TCGA-A8-A09I-01'

    print('load counts and PSI')

    counts = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/counts.TCGA-BRCA.txt',sep='\t')

    events = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-filtered-names-75p.txt',sep='\t')

    print('start to run')

    tumorPSI = events[id_].values

    tumorPSI[np.isnan(tumorPSI)] = 0.0

    tumorPSI = tumorPSI.tolist()


    UID = events['UID'].tolist()

    part = []
    for i in range(len(UID)):
        uid = UID[i]
        x = uid.split('|')
        try: x[0].split(':')[3]
        except IndexError: event = x[0].split(':')[2]
        else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])
        query = uid.split('|')[1].split(':')[0] + ':' + event
        part.append(query)
        
    counts_id = counts['AltAnalyze_ID'].tolist()

    counts_id_useful = [item.split('=')[0] for item in counts_id]


    new = [name[:15] for name in counts.columns]
    counts.columns = new


    query_id_counts = counts[id_].tolist()



        
    dic = {}
    for i in range(len(counts_id_useful)):
        ensg = counts_id_useful[i].split(':')[0]
        event = counts_id_useful[i].split(':')[1:]   
        event = ''.join(event)
        info = query_id_counts[i]
        try:
            dic[ensg].append((event,info))
        except KeyError:
            dic[ensg] = []
            dic[ensg].append((event,info))



    counts_extract = []
    for i in range(len(part)):
        flag = False
        part_ensg = part[i].split(':')[0]
        searchSpace = dic[part_ensg]
        for j in searchSpace:
            match = part_ensg+':'+j[0]
            if part[i] == match:
                counts_extract.append(j[1])
                flag = True
                break
        if flag == False:
            counts_extract.append(0)


    final = pd.DataFrame({'UID':UID,'PSI':tumorPSI,'count':counts_extract})

    cond = []
    for i in range(final.shape[0]):
        psi = final['PSI'].iloc[i]
        if psi == 0.0:
            cond.append(False)
        else:
            cond.append(True)
    final['cond'] = cond
    final = final[final['cond']]
    final = final.drop(columns=['cond'])
    final = final.set_index(pd.Index(np.arange(final.shape[0])))

    final.to_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/{0}.ori.txt'.format(id_),sep='\t',index=None)