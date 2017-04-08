
#command: python defects_production.py vcluster_production_reactor.txt


import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=24)


class defects_production:
    def __init__(self,filename,scaling_factor=1.0):
        self.cluster_spectrum = np.genfromtxt(filename,skip_header=1)
        self.cluster_spectrum[:,1] *= scaling_factor

    def printInfo(self):
        print 'size   #/um^3s'
        first = []
        second = []
        total = 0.0
        for row in self.cluster_spectrum:
            if row[1] > 1.0e-10:
                first.append(int(row[0]))
                second.append(row[1])
                total += first[-1]*second[-1]
        print str(first).replace(",","")
        #print str(['%.4f' %x for x in second ]).replace("'","").replace(",","")
        print str(['%d' %int(x) for x in second ]).replace("'","").replace(",","")
        print 'Total point defects: ', total




if __name__ == '__main__':
    #filename = 'vcluster_production_reactor.txt'
    filename = sys.argv[1]

    #1. neutron irradiation at 90C, 2.28e-7 dpa/s
    #ref: Defect evolution in single crystalline tungsten following low temperature and low dose neutron irradiation
    #scaling = 2.28e-7/1.5825e-11 #dpa/s / atomic_volume_um^3


    #2. ion irradiation at 30C
    #ref: Direct oberservation of size scaling 
    scaling = 7.9e8 #150keV 0.0125dpa/s
    #scaling = 1.02e9 #400keV 0.01613dpa/s
    

    print "Read data from " , filename
    d = defects_production(filename,scaling_factor=scaling)
    d.printInfo()
