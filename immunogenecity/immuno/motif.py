import os
import pandas as pd

def hla_convert(hla):
    hla = hla.replace('*','')
    return hla

def write_fasta(pep,hla,label):
    hla = hla_convert(hla)
    if not os.path.exists('motif/{0}'.format(hla)): os.makedirs('motif/{0}'.format(hla))
    with open('motif/{0}/{0}_{1}_9.fasta'.format(hla,label),'w') as f1,open('motif/{0}/{0}_{1}_10.fasta'.format(hla,label),'w') as f2:
        for p in pep:
            if len(p) == 9:
                f1.write('>placeholder\n')
                f1.write('{0}\n'.format(p))
            elif len(p) == 10:
                f2.write('>placeholder\n')
                f2.write('{0}\n'.format(p))


if __name__ == '__main__':

ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')
hla_score = ['HLA-A*0203','HLA-A*0205','HLA-A*0206','HLA-A*0207','HLA-A*3001','HLA-A*6801','HLA-A*9235','HLA-A*9253','HLA-B*0602','HLA-B*1557','HLA-B*1801','HLA-B*3505','HLA-B*4001','HLA-B*4002','HLA-B*4044','HLA-B*4102','HLA-B*4104','HLA-B*4202','HLA-B*4601','HLA-B*5706','HLA-B*5712','HLA-B*8103','HLA-B*9234','HLA-C*0102','HLA-C*0304','HLA-C*0401','HLA-C*0517','HLA-C*0602','HLA-C*0756','HLA-C*1604','HLA-A*0101','HLA-A*0201','HLA-A*0224','HLA-A*0301','HLA-A*0362','HLA-A*1101','HLA-A*2301','HLA-A*2402','HLA-A*3003','HLA-A*9234','HLA-B*0702','HLA-B*0801','HLA-B*1402','HLA-B*1501','HLA-B*2703','HLA-B*2705','HLA-B*2709','HLA-B*2713','HLA-B*3501','HLA-B*3508','HLA-B*3901','HLA-B*4201','HLA-B*4402','HLA-B*4403','HLA-B*4405','HLA-B*5101','HLA-B*5301','HLA-B*5701','HLA-B*5703','HLA-B*5801','HLA-B*8102','HLA-C*1510']

tmp1,tmp2 = ori.groupby(by=['immunogenecity'])
negative = tmp1[1]
positive = tmp2[1]

for hla in hla_score:
    pep = negative.loc[negative['HLA']==hla]['peptide']
    write_fasta(pep,hla,'negative')

for hla in hla_score:
    pep = positive.loc[positive['HLA']==hla]['peptide']
    write_fasta(pep,hla,'positive')



