import numpy as np
import scipy as sp
from scipy import linalg as lig
from scipy import stats as stats

def Mod_matrix(A):  #calculate mod of each vector in a metrix
    B = np.array(A)
    for i in range(A.shape[1]):
        B[:,i] = A[:,i]/lig.norm(A[:,i])
    return B

def Minus_mod(A):  # generate n-1 vectors from n sample points
    E = np.array(A[:,:-1])
    for i in range(1,A.shape[1]):
        E[:,i-1] = A[:,i] - A[:,i-1]
    return Mod_matrix(E)
    
def Angle_sector(E): # calculate the angle secoter of a series vectors
    A = np.dot(np.logspace(0,0,E.shape[1]), lig.pinv((np.dot(sp.transpose(E),E))))
    B = np.dot(E,sp.transpose(A))
    B = B/lig.norm(B)
    return B

def Best_direction(A): # calculate the co-bisectors
    E = Minus_mod(A)
    return Angle_sector(E)

def Angle(A): # calculate the angel berween bisector and vector
    E = Best_direction(A)
    F = Minus_mod(A)
    return np.arccos(np.dot(sp.transpose(F),E))

def Project(A,e): # projecte sample matrix A to a vector e
    return np.dot(sp.transpose(e),A)

def Read_Matrixformat(Filename): # Read GEO SeriesMatrix formate data, generate a line-list
    File_li = open(Filename).readlines()
    File_list = []
    for line in File_li:
        if '!series_matrix_table_begin' in line:
            File_list = File_li[File_li.index(line)+1:-1]
            break
    return File_list
            
def File2narray(File_list,x): # transform line-list to a narray, x is 1 when only Probe_ID, 
                                                      # data will be sorted and meaned acoording to Order_file,
                                                      # generate 1 narray: Ordered_data and 1 list: Probe_list
    print 'Transforming data to Narray object...\n'
    Sample = File_list[0][:-1].split('\t')[x:]  #generate sample list
    print Sample
    Linenu = len(File_list)-1  #Linenu is the number of Probe    
    Colnu = len(Sample)  #Colnu is the number of samples
    print 'The file contain ', Colnu, ' samples;', Linenu,' Probe\n'
    print 'Data transforming...\n'
    Chipdata = np.zeros(Linenu*Colnu, dtype= 'f').reshape(Linenu,Colnu)  #generate probe*sample narray ,
                                                                                                         # value is 1
    Probe_list = []
    for i in range(1,Linenu+1):   #denote each expression index to right position in narray object
        Line = File_list[i]
        for j in range(x,Colnu+x):
            #print i, j 
            Element = Line.split('\t')
            Chipdata[i-1,j-x] = Element[j]
        Probe_list.append(Element[0])
    return Chipdata, Probe_list, Sample

def Isolate_order(Chipdata,Order_file,Sample):
    print 'Sorting data accorfing to Order_file...\n'
    Linenu = Chipdata.shape[0]
    Order_list = open(Order_file).readlines()
    Statenu = len(Order_list)
    Ordered_data = np.logspace(0,0,Linenu*Statenu).reshape(Linenu,Statenu)
    n = 0
    for i in range(Statenu):
        name = Order_list[i][(Order_list[i].index(':'))+1:-1]
        Rep = name.split(',')
        print 'State ',i,' have ', len(Rep),' Samples: '
        for j in range(len(Rep)):
            print Rep[j]
            Ordered_data[:,i] = Ordered_data[:,i] + Chipdata[:,int(Sample.index(Rep[j]))]
        Ordered_data[:,i] = Ordered_data[:,i]/(j+1)
    return Ordered_data 
    

def Median_normalize(data):
    for i in range(data.shape[1]):
        data[:,i] = data[:,i] - np.median(data[:,i])
    return data

def Log2_data(data):
    return np.log2(data)

def Data_statis(data):
    Colnu = data.shape[1]
    Data_statis = np.logspace(0,0,Colnu*4).reshape(4,Colnu)  #Make a narray to store the max,
                                                                                        #mean, min, std of Chipdata 
    for i in range(Colnu):
        Data_statis[0,i] = data[:,i].max()
        Data_statis[1,i] = data[:,i].min()
        Data_statis[2,i] = data[:,i].mean()
        Data_statis[3,i] = data[:,i].std()
    print 'The max in data: ',  Data_statis[0,:]
    print 'The min in data: ',  Data_statis[1,:]
    print 'The mean in data: ',  Data_statis[2,:]
    print 'The std in data: ',  Data_statis[3,:]

def P_FDR(Weight,Probe_list, FDR):  # Weight come from best_direction, Probe_list come form
                                                       # File2narray. FDR is input. Calculate P-value and modified P-value, 
                                                       # write Probe, weight, P-value to file Significant candi.txt
    Table_type = np.dtype([('Probe',np.str_,16),('Weight',float),('P-value',float),('FDR',float)])
    P_table = np.zeros(len(Probe_list),dtype=Table_type)
    print 'Calculating P-value...\n'
    Pdf = stats.norm.pdf(Weight,np.mean(Weight),np.std(Weight))
    Pdf = Pdf*0.01
    print np.shape(P_table)
    Fdic = {}
    for i in range(len(Probe_list)):
        P_table[i]['Probe'] = Probe_list[i]
        P_table[i]['Weight'] = Weight[i]
        P_table[i]['P-value'] = Pdf[i]
        Fdic[Weight[i]*len(Probe_list)+(i/float(len(Probe_list))/len(Probe_list))] = Probe_list[i]+'\t'+ str(Weight[i])+'\t'+str(Pdf[i])
    Sigout = file('Significant candi.txt','w')
    Sigout.write('Probe     Weight  P   BH-P(FDR< '+str(FDR)+')\n')
    x = 1
    r = int(0.5*len(Probe_list)-2)
    print len(Probe_list),len(Fdic)
    while x < r:
        weight = max(Fdic)
        Probe_Pvalue = Fdic.pop(weight)
        pvalue = float(Probe_Pvalue.split('\t')[-1])
        K_hat = 2*x*FDR/(len(Probe_list))
        #if K_hat > pvalue:
        BHP = pvalue*(len(Probe_list)/x)
        Sigout.write(Probe_Pvalue+'\t'+str(BHP)+'\n')
        x = x+1

    x = 1
    while x < r:
        weight = min(Fdic)
        Probe_Pvalue = Fdic.pop(weight)
        pvalue = float(Probe_Pvalue.split('\t')[-1])
        K_hat = 2*x*FDR/(len(Probe_list))
        #if K_hat > pvalue:
        BHP = pvalue*(len(Probe_list)/x)
        Sigout.write(Probe_Pvalue+'\t'+str(BHP)+'\n')
        x = x+1
        
    Sigout.close()





f1=open('GSE21630_data').readlines()
(chipdata,Pr_list,sample) = File2narray(f1,1)

data1 = Isolate_order(chipdata,'B-cell-mir-order',sample)


data1 = Log2_data(data1)


data1 = Median_normalize(data1)

ou = file('refined_exp','w')
for i in range(data1.shape[0]):
    ou.write(Pr_list[i]+'\t'+str(data1[i,:])+'\n')
ou.close()

Data_statis(data1)


Weight1 = Best_direction(data1)
Ptable1 = P_FDR(Weight1, Pr_list, 1)

m1 = Project(data1,Weight1)
print 'The sample 1 on 1st line project : ',m1-m1[0]


print Angle(data1), Angle(data1)*180/np.pi
