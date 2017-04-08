import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=20)
from scipy.optimize import curve_fit


class defects:
    #read data, e.g. energy intervals from file
    def __init__(self,filename):
        self.pka_spectrum = np.genfromtxt(filename,skip_header=1)
        self.pka_spectrum = self.pka_spectrum[self.pka_spectrum[:,2].argsort()]
        self.pka_count = np.zeros(len(self.pka_spectrum)) #count the number of points in each interval

    #get the proportional value for this energy
    def getProp(self,energy):
        if energy >= self.pka_spectrum[-1,1]:
            return self.pka_spectrum[-1,-1]
        for line in self.pka_spectrum:
            if energy>=line[0] and energy<line[1]:
                return line[-1]
        return 0.0

    #set the number of pkas in each energy interval
    def setCount(self,energy):
        if energy >= self.pka_spectrum[-1,1]:
            self.pka_count[-1] += 1
        for i,line in enumerate(self.pka_spectrum):
            if energy>=line[0] and energy<line[1]:
                self.pka_count[i] += 1
                break

    #get the number of pkas in each energy interval
    def getCount(self,energy):
        if energy >= self.pka_spectrum[-1,1]:
            return self.pka_count[-1]
        for i,line in enumerate(self.pka_spectrum):
            if energy>=line[0] and energy<line[1]:
                return self.pka_count[i]

    #plot number of defects vs. pka energy, and defect cluster distribution
    def plot_defects(self,filenames,energies,plotPointDefects=False,plotDist=False):
        data = dict() # energy: cluster distribution
        total_pointdefects = []
        for i,filename in enumerate(filenames):
            #check existence of file
            if os.path.exists(filename) and os.path.getsize(filename) > 0:
                tmp = np.genfromtxt(filename).reshape(-1,2)
                data[str(energies[i])] = tmp
                num = np.sum(tmp[:,0]*tmp[:,1]) 
                total_pointdefects.append(num)

        if plotPointDefects:
            plt.subplot(111)
            plt.plot(energies,total_pointdefects,'ro-',markersize=5,linewidth=4,label='MD annealing 30K')
            plt.xscale('log')
            plt.xlabel('PKA energy (keV)')
            plt.ylabel('Number of point defects')
            #plt.legend(loc='upper left', shadow=True,prop={'size':14})
            plt.tight_layout()
            plt.grid()
            plt.savefig('pka_pointdefects.png',dpi=300)
            plt.close()
            
        if plotDist: 
            plt.subplot(111)
            accumulate = dict()
            for e in energies:
                self.setCount(e) #set count of energies in each interval
            for e in energies:
                prop = self.getProp(e)
                count = self.getCount(e)
                for line in data[str(e)]:
                    key = line[0]
                    if key in accumulate:
                        accumulate[key] += line[1]*prop/count
                    else:
                        accumulate[key] = line[1]*prop/count
            max_size = np.max(np.array(accumulate.keys()))
            cluster_data = []
            avg_cluster_size = 0.0 #mean
            for key in range(1,int(max_size)+1):
                if key in accumulate:
                    cluster_data.append([accumulate[key]])
                else:
                    cluster_data.append([0.0])
                avg_cluster_size += key*1.0*cluster_data[-1][0]

            f = open('icluster_production_400keV.txt','w')
            f.write('Size   Number\n')
            for i,item in enumerate(cluster_data):
                print item[0], avg_cluster_size
                print >> f,'%d' %(i+1),'%.8f' %(item[0]/avg_cluster_size)
            f.close()
    
            dim = len(cluster_data[0])
            dimw = 1 #width of a pillar
            w = dimw * dim #width of a field
            x = np.arange(len(cluster_data))
            colors = matplotlib.cm.jet(np.linspace(0,1,dim))
            for i in range(dim) :
                y = [d[i] for d in cluster_data]
                b = plt.bar(x + i * dimw, y, dimw, bottom=0.001,label='20 random directions for each PKA energy',color='r')
            num_ticks = 10
            int_ticks = len(x)/10 if len(x)/10>=1 else 1
            ticks = list(x[0:-1:int_ticks]+0.5)
            plt.xticks(ticks, fontsize = 16, fontweight='bold')
            labels = [str(x[i]+1) for i in range(0,len(x),int_ticks)]
            plt.gca().set_xticklabels(labels)
            plt.ylabel('Number of cluster')
            plt.xlabel('Cluster size')
            plt.legend(loc='upper right', shadow=True,prop={'size':14,'weight':'bold'})
            plt.tight_layout()
            plt.savefig('icluster_hist_400keV.png',dpi=300)
            plt.close()
        
                
def simplePlot():
    data_30K = np.genfromtxt('30K_Epka_Ndefects.txt',delimiter=',')
    #data_363K = np.genfromtxt('363K.txt',delimiter=',')
    plt.subplot(111)
    plt.plot(data_30K[0,:],data_30K[1,:],'ro-',markersize=5,linewidth=4,label='MD annealing 30K')
    #plt.plot(data_363K[0,:],data_363K[1,:],'bo-',markersize=5,linewidth=4,label='MD annealing 363K')
    plt.xscale('log')
    plt.xlabel('PKA energy (keV)')
    plt.ylabel('Number of point defects')
    plt.legend(loc='upper left', shadow=True,prop={'size':14,'weight':'bold'})#
    plt.tight_layout()
    plt.grid()
    plt.savefig('pka_pointdefects.png',dpi=300)
    plt.close()



if __name__ == '__main__':
    '''
    filename = 'pka_400keV.txt'
    #filename = 'pka_400keV.txt'
    #energies = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1,  1.5, 2,  3,  5,  7.5, 10, 15, 20, 30, 40, 50, 60, 75, 100, 120, 150]
    energies = [0.1,0.16,0.25,0.4,0.63,1, 1.58,2.51,4, 6.3,10,15.8, 25.1,40,63, 100, 150, 250]
    d = defects(filename)
    filenames = []
    for e in energies:
        filenames.append(str(e)+'keV/icluster_frequency.txt')
    d.plot_defects(filenames,energies,plotPointDefects=True,plotDist=True)
    '''
    simplePlot()
