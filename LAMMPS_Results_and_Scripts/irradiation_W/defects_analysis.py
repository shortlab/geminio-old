import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=20)
from scipy.optimize import curve_fit


class IrradiationPKA:
    #read data, e.g. energy intervals from file : lower limit (keV); higher limit (keV); mean (keV); proportion; variance
    def __init__(self,data):
        self.pka_spectrum = np.copy(data)
            

class defects:
    def __init__(self,pka):
        self.pka = pka
        self.pka_count = np.zeros(len(self.pka.pka_spectrum)) #count the number of points in each interval

    #get the proportional value for this energy
    def getProp(self,energy):
        if energy >= self.pka.pka_spectrum[-1,1]:
            return self.pka.pka_spectrum[-1,3]
        for line in self.pka.pka_spectrum:
            if energy>=line[0] and energy<line[1]:
                return (line[3],line[4]) #prop, variance
        return (0.0,0.0)

    #set the number of pkas in each energy interval
    def setCount(self,energy):
        if energy >= self.pka.pka_spectrum[-1,1]:
            self.pka_count[-1] += 1
        for i,line in enumerate(self.pka.pka_spectrum):
            if energy>=line[0] and energy<line[1]:
                self.pka_count[i] += 1
                break

    #get the number of pkas in each energy interval
    def getCount(self,energy):
        if energy >= self.pka.pka_spectrum[-1,1]:
            return self.pka_count[-1]
        for i,line in enumerate(self.pka.pka_spectrum):
            if energy>=line[0] and energy<line[1]:
                return self.pka_count[i]

    #plot number of defects vs. pka energy, and defect cluster distribution
    def plot_defects(self,filenames,energies,plotPointDefects=False,plotDist=False,writeToFile=None):
        data = dict() # energy: cluster distribution
        total_pointdefects = []
        for i,filename in enumerate(filenames):
            #check existence of file
            if not os.path.exists(filename) :
                print 'filename : ', filename
                raise Exception('File not exist')
            if os.path.getsize(filename) == 0 :
                tmp = np.array([[1,0]]) #empty file
            else :
                ###!tmp = np.genfromtxt(filename,skip_header=1).reshape(-1,2) #Size Number Variance
                tmp = np.genfromtxt(filename).reshape(-1,2) #Size Number
            data[str(i)] = tmp
            num = np.sum(tmp[:,0]*tmp[:,1]) 
            total_pointdefects.append(num)

        if plotPointDefects:
            plt.subplot(111)
            plt.plot(energies,total_pointdefects,'ro-',markersize=5,linewidth=4,label='MD annealing 30K')
            plt.xscale('log')
            plt.xlabel('PKA energy (keV)')
            plt.ylabel('Number of point defects')
            plt.legend(loc='best', shadow=True,prop={'size':14,'weight':'bold'})
            plt.tight_layout()
            plt.grid()
            plt.savefig('pka_pointdefects.png',dpi=300)
            plt.close()
            
        if plotDist: 
            plt.subplot(111)
            accumulate = dict()
            ###!avg_squared = dict()
            for e in energies:
                self.setCount(e) #set count of energies in each interval
            for i,e in enumerate(energies):
                (prop,variance) = self.getProp(e)
                count = self.getCount(e)
                for line in data[str(i)]:
                    key = int(line[0])
                    if key in accumulate:
                        accumulate[key] += line[1]*prop/count
                        ###!avg_squared[key] += ((variance+prop**2)*(line[2]+line[1]**2)-(line[1]*prop)**2)/count**2
                    else:
                        accumulate[key] = line[1]*prop/count
                        ###!avg_squared[key] = ((variance+prop**2)*(line[2]+line[1]**2)-(line[1]*prop)**2)/count**2
            max_size = np.max(np.array(accumulate.keys()))
            cluster_data = []
            ###!var = []
            avg_cluster_size = 0.0 #mean
            for key in range(1,int(max_size)+1):
                if key in accumulate:
                    cluster_data.append([accumulate[key]])
                    ###!var.append(avg_squared[key])
                else:
                    cluster_data.append([0.0])
                    ###!var.append(0.0)
                avg_cluster_size += key*1.0*cluster_data[-1][0] #\sum_i i*f(i) total number of point defects

            #trim data from where variance is very large
            end_index = len(cluster_data)
            ###!for i,item in enumerate(cluster_data):
            ###!    if item[0]>1.0e-10 : 
            ###!        if var[i]>0 and float(item[0])/np.sqrt(var[i])<1.4 : 
            ###!            end_index = i
            ###!            break
            ###!for i in range(end_index,len(cluster_data)):
            ###!    avg_cluster_size -= (i+1)*cluster_data[i][0]
                        
            f = open(writeToFile+'.txt','w')
            f.write('Size   Number   Std\n')
            for i in range(end_index):
                if cluster_data[i][0]>1.0e-10 : 
                    ###!print >> f,'%d' %(i+1),'%.8f' %(cluster_data[i][0]/avg_cluster_size), '%.8f' %(np.sqrt(var[i])/avg_cluster_size)
                    print >> f,'%d' %(i+1),'%.8f' %(cluster_data[i][0]/avg_cluster_size)
            f.close()
    
            dim = len(cluster_data[0])
            dimw = 1 #width of a pillar
            w = dimw * dim #width of a field
            x = np.arange(len(cluster_data))
            colors = matplotlib.cm.jet(np.linspace(0,1,dim))
            for i in range(dim) :
                y = [d[i] for d in cluster_data]
                b = plt.bar(x + i * dimw, y, dimw, bottom=0.001,label='30 random directions for each PKA energy',color='b')
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
            plt.savefig(writeToFile+'.png',dpi=300)#'icluster_hist_400keV.png'
            plt.close()
        
                    



if __name__ == '__main__':
    
    part1 = 1
    part2 = 0
    if part1 :
        #part 1: ion irradiation
        filename = 'pka_400keV_1e5ions.txt'
        data = np.genfromtxt(filename,skip_header=1)
        pka_spectrum = data
        mean_energies = list(data[:,2]) #real pka energy
        MD_mean_energies = [0.1 ,0.14 ,0.22 ,0.35 ,0.55 ,0.86 ,1.34 ,2.1 ,3.3 ,5.13 ,8.02 ,12.5 ,19.5 ,30 ,46.9 ,72.7 ,107 ,171.9] #pka energy partitioned to displace atoms
        ifilenames = []
        for e in MD_mean_energies:
            ifilenames.append('30K'+str(e)+'keV/icluster_frequency.txt')
        vfilenames = []
        for e in MD_mean_energies:
            vfilenames.append('30K'+str(e)+'keV/vcluster_frequency.txt')

    if part2 :
        #part 2: neutron irradiation
        filename = 'pka_reactor.csv'
        data = np.genfromtxt(filename,skip_header=1,delimiter=',')
        low = []
        high = []
        mean = []
        prop = []
        for i in range(11,len(data)):
            low.append(data[i,2])
            high.append(data[i,3])
            mean.append(data[i,4])
            prop.append(data[i,6])
        pka_spectrum = np.transpose(np.array([low,high,mean,prop]))
        mean_energies = mean #real pka energies
        MD_mean_energies=[0.1 ,0.15 ,0.2 ,0.35 ,0.63 ,0.75 ,1.34 ,2.1 ,3 ,4 ,6.3 ,8.02 ,10 ,12.6 ,19.5 ,23 ,30.4 ,45 ,55 ,72.7 ,100] #pka energy partitioned to displace atoms
        ifilenames = []
        for e in MD_mean_energies:
            ifilenames.append(str(e)+'keV/icluster_frequency.txt')
        vfilenames = []
        for e in MD_mean_energies:
            vfilenames.append(str(e)+'keV/vcluster_frequency.txt')
        filenames = []

    if len(mean_energies) != len(MD_mean_energies) :
        raise Exception('Number of data does not match')

    PKA = IrradiationPKA(pka_spectrum)
    d = defects(PKA)
    d.plot_defects(ifilenames,mean_energies,plotPointDefects=True ,plotDist=True,writeToFile='icluster_production_400keV')
    d.plot_defects(vfilenames,mean_energies,plotPointDefects=False,plotDist=True,writeToFile='vcluster_production_400keV')
