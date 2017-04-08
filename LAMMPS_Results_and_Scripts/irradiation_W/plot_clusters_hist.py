import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=20)


def plotHist(data,labels,tick_labels,outfilename='results'):
    #plot histogram for energy spectrum
    plt.subplot(111)
    dim = len(data[0])
    dimw = 0.5 #width of a pillar
    w = dimw * dim #width of a field
    x = np.arange(len(data))
    colors = matplotlib.cm.gnuplot(np.linspace(0,1,dim))
    for i in range(dim) :
        y = [d[i] for d in data]
        b = plt.bar(x + i * dimw, y, dimw, bottom=0.001,label=labels[i],color=colors[i])
    #ticks = list(x[0:-1:2]+w/2)
    #plt.xticks(ticks, fontsize = 16, fontweight='bold')
    plt.gca().set_xticks(x+w/2)
    plt.gca().set_xticklabels(tick_labels)

    plt.ylabel('Proportion')
    plt.xlabel('Cluster size')
    plt.legend(loc='upper right', shadow=True,prop={'size':14,'weight':'bold'})
    plt.tight_layout()
    plt.savefig(outfilename+'.png',dpi=300)

def plotClusterProduction():
    filenames = ['icluster_production_150keV.txt','vcluster_production_150keV.txt','icluster_production_400keV.txt','vcluster_production_400keV.txt']
    markers = ['o','s','o','s']
    colors = ['r','r','b','b']
    markerfacecolors = ['w','w','b','b']
    edgecolors = ['r','r','b','b']
    labels = ['SIA-150keV','Vacancy-150keV','SIA-400keV','Vacancy-400keV']
    plt.subplot(111)
    for i,filename in enumerate(filenames) :
        data = np.genfromtxt(filename,skip_header=1)
        #value_error_ratio = list(data[:,1]/data[:,2]<1.1)
        #end_index = next((i for i, v in enumerate(value_error_ratio) if v != False), len(value_error_ratio))
        end_index = len(data)
        plt.plot(data[:end_index,0],data[:end_index,1],marker=markers[i],markersize=5,mew=2,markeredgecolor=edgecolors[i],markerfacecolor=markerfacecolors[i],linestyle='None',label=labels[i])#,color=colors[i]
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Proportion')
    plt.xlabel('Cluster size')
    plt.legend(loc='best', shadow=True,prop={'size':14,'weight':'bold'})#
    plt.grid()
    plt.tight_layout()
    plt.savefig('cluster_production.png',dpi=300)


    pass

if __name__ == '__main__':
    plot_hist = False
    if plot_hist :
        filename1 = 'icluster_production_150keV.txt'
        idata_150keV = np.genfromtxt(filename1,skip_header=1)
        filename2 = 'vcluster_production_150keV.txt'
        vdata_150keV = np.genfromtxt(filename2,skip_header=1)
        
        isizes = list(idata_150keV[:,0].astype('int'))
        vsizes = list(vdata_150keV[:,0].astype('int'))
        allsizes = list(isizes)
        [allsizes.append(i) for i in vsizes if i not in isizes]
        
        allsizes = sorted(allsizes,key=int)
       
        idata = []
        vdata = []
        for i in allsizes :
            if i in isizes :
                idata.append(idata_150keV[isizes.index(i),1])
            else :
                idata.append(0)
            if i in vsizes :
                vdata.append(vdata_150keV[vsizes.index(i),1])
            else :
                vdata.append(0)
        data = [ [idata[i],vdata[i]] for i in range(len(allsizes))]
    
        labels = ['Interstitial','Vacancy']
        tick_labels = [str(i) for i in allsizes]
        plotHist(data,labels,tick_labels,outfilename='cluster_production_150keV')

    plot_lines = True
    if plot_lines :
        plotClusterProduction()
