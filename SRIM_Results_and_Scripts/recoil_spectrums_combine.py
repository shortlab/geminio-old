import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=20)
from scipy.optimize import curve_fit

def pka_spectrum(filename,output_name='result'):
    
    #define energy bins in keV
    start = 2
    end = 6
    step = 0.2
    exponents = np.arange(start,end,step)*1.0
    rho = np.zeros(len(exponents)) #histogram
    rho_sq = np.zeros(len(exponents)) #histogram for variance
    total_pkas = 0.0
    left_ends = []
    right_ends = []
    medians = []
    for i in range(len(exponents)):
        right = 10**(exponents[i]+step/2.0)
        left = 10**(exponents[i]-step/2.0)
        left_ends.append(left/1000.0)
        right_ends.append(right/1000.0)
        medians.append(10**(exponents[i]*1.0)/1000.0)


    #recoil spectrum for every ions
    total_pkas = 0 #keep track of the total pkas created
    total_ions = 0
    with open(filename) as f:
        while 1 :
            line = f.readline()
            if not line:
                break
            if line[0:5] == '-----' :
                pka_energy = []
                while 1 : 
                    line = f.readline()
                    if not line:
                        break
                    if line[0:5] != '=====':
                        pka_energy.append(float(line.split()[7]))
                    else :
                        break
                if len(pka_energy) == 0:
                    continue
                total_pkas += len(pka_energy)
                total_ions += 1
                pka_energy = np.array(pka_energy)
                for i in range(len(exponents)):
                    right = 10**(exponents[i]+step/2.0)
                    left = 10**(exponents[i]-step/2.0)
                    current_count = np.sum((pka_energy<=right) & (pka_energy>left))
                    rho[i] += current_count

    rho  = rho/total_pkas #proportion
    rho_sq = rho*(1.0-rho)/total_pkas #sample proportion variance  
    print 10**exponents/1000.0
    print 'Proportion: '
    print rho
    print 'total pka: ', total_pkas
    print 'total ions: ', total_ions
    print 'total proportion (correct if 1.0): ', np.sum(rho)
    print 'Write to file: ', output_name+'.txt'
    f = open(output_name+'.txt','w')
    f.write('Left[keV]  right[keV]  median[keV]  Prop  Variance\n')
    for i in range(len(exponents)):
        print >> f,'%.4f' %left_ends[i],'%.4f' %right_ends[i],'%.4f' %medians[i],rho[i],rho_sq[i]
    f.close()

    return (exponents,list(rho))


def plotHist(data,labels,tick_labels,filename='results'):
    #plot histogram for energy spectrum
    plt.subplot(111)
    dim = len(data[0])
    dimw = 0.4 #width of a pillar
    w = dimw * dim #width of a field
    x = np.arange(len(data))
    colors = matplotlib.cm.gnuplot(np.linspace(0,1,dim))
    for i in range(dim) :
        y = [d[i] for d in data]
        b = plt.bar(x + i * dimw, y, dimw, bottom=0.0,label=labels[i],color=colors[i])
    ticks = list(x[0:-1:5]+w/2)
    print x
    print ticks
    plt.xticks(ticks, fontsize = 20)
    #plt.gca().set_xticks(x+w/2)
    plt.gca().set_xticklabels(tick_labels)

    plt.ylabel('Proportion')
    plt.xlabel('PKA energy (keV)')
    plt.legend(loc='upper right', shadow=True,prop={'size':14,'weight':'bold'})
    plt.tight_layout()
    plt.savefig(filename+'_hist.png',dpi=300)
    plt.close()
    
                    



if __name__ == '__main__':
    '''
    filename1 = 'SRIM_Outputs_150keV_1e4ions/COLLISON.txt'
    exponents,data_150keV = pka_spectrum(filename1,output_name='pka_150keV_1e4ions')
    filename2 = 'SRIM_Outputs_400keV_1e4ions/COLLISON.txt'
    exponents,data_400keV = pka_spectrum(filename2,output_name='pka_400keV_1e4ions')

    labels = ['150keV','400keV']
    tick_labels = [str('$\mathbf{10^{'+'%g' %(exponents[i]-3)+'}}$') for i in range(0,len(exponents),2)]
    data  = [ [data_150keV[i], data_400keV[i] ] for i in range(len(data_150keV))]
    plotHist(data,labels,tick_labels)
    #labels = '10000 150 keV W statistics'
    '''

    filename1 = 'SRIM_Outputs_150keV_1e5ions/COLLISON.txt'
    #filename1 = 'SRIM_Outputs_150keV/COLLISON.txt'
    exponents,data_150keV = pka_spectrum(filename1,output_name='pka_150keV_1e5ions')
    filename2 = 'SRIM_Outputs_400keV_1e5ions/COLLISON.txt'
    #filename2 = 'SRIM_Outputs_400keV/COLLISON.txt'
    exponents,data_400keV = pka_spectrum(filename2,output_name='pka_400keV_1e5ions')

    labels = ['150keV','400keV']
    #tick_labels = [str('$\mathbf{10^{'+'%g' %(exponents[i]-3)+'}}$') for i in range(0,len(exponents),2)]
    tick_labels = ['$\mathbf{10^{'+str('%d' %(i))+'}}$' for i in range(-1,3)]
    data  = [ [data_150keV[i], data_400keV[i] ] for i in range(len(data_150keV))]
    plotHist(data,labels,tick_labels)
