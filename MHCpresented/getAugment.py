import pandas as pd
import ast
import pickle
neo_antigen8 = pd.read_csv('NeoJunction_8_mark.txt',sep='\t',usecols=['MHCIresult'])
neo_antigen_list8 = neo_antigen8['MHCIresult'].tolist()
neo_antigen9 = pd.read_csv('NeoJunction_9_mark.txt',sep='\t',usecols=['MHCIresult'])
neo_antigen_list9 = neo_antigen9['MHCIresult'].tolist()
neo_antigen10 = pd.read_csv('NeoJunction_10_mark.txt',sep='\t',usecols=['MHCIresult'])
neo_antigen_list10 = neo_antigen10['MHCIresult'].tolist()
neo_antigen11 = pd.read_csv('NeoJunction_11_mark.txt',sep='\t',usecols=['MHCIresult'])
neo_antigen_list11 = neo_antigen11['MHCIresult'].tolist()

def getb(list_):
    sb,wb = [],[]
    for item in list_:
        if item == 'No candidates': continue
        else: 
            info = ast.literal_eval(item)
            pep = info['HLA-A29:02']
            pep_sb,pep_wb = pep[0],pep[1]
            sb.extend(pep_sb)
            wb.extend(pep_wb)
    b = sb + wb
    return b
   
b8 = getb(neo_antigen_list8)
b9 = getb(neo_antigen_list9)
b10 = getb(neo_antigen_list10)
b11 = getb(neo_antigen_list11)
btotal = b8+b9+b10+b11
len(btotal)
type(btotal)

# write to fasta
def toFasta(list_):
    with open('augment.fasta','w') as file1:
        for index,item in enumerate(list_):
            file1.write('>mer{0}\n'.format(index+1))
            file1.write('{0}\n'.format(item))
            
toFasta(b)
