from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
import numpy as np
import sys
import itertools
import collections

#command: ./ovitos ovito_process.py dump_irra_minimized.data dump_minimized.100284.data
#sys.argv[1]: data file to analyze
#sys.argv[2]: reference file to perform wigner-seitz analysis

def cluster_analysis(filename,reference_filename,boxL):
    node = import_file(filename)
    
    # Perform Wigner-Seitz analysis:
    ws = WignerSeitzAnalysisModifier(
        per_type_occupancies = True, 
        eliminate_cell_deformation = True)
    ws.reference.load(reference_filename)
    node.modifiers.append(ws)
    node.compute()
    print("Number of vacant sites: %i" % node.output.attributes['WignerSeitz.vacancy_count'])
    print("Number of interstitial: %i" % node.output.attributes['WignerSeitz.interstitial_count'])
    # Define a modifier function that selects sites of type A=1 which
    # are occupied by exactly one atom of type B=2.
    def modify(frame, input, output):
    
        # Retrieve the two-dimensional Numpy array with the site occupancy numbers.
        occupancies = input.particle_properties['Occupancy'].array
        
        # Get the site types as additional input:
        site_type = input.particle_properties.particle_type.array
    
        # Set up a particle selection by creating the Selection property:
        selection = output.create_particle_property(ParticleProperty.Type.Selection).marray
        
        # Note that the Occupancy array uses 0-based
        # indexing, while atom type IDs are typically 1-based.
        selection[:] = (site_type == 1) & (occupancies[:] < 1)
        
        # Additionally output the total number of defects as a global attribute:
        output.attributes['Defect_count'] = np.count_nonzero(selection)
    
    # Insert Python modifier into the data pipeline.
    node.modifiers.append(PythonScriptModifier(function = modify))
    
    node.modifiers.append(InvertSelectionModifier())
    node.modifiers.append(DeleteSelectedParticlesModifier())
    node.compute()
    
    lat_param = 3.1658
    coords = node.output.particle_properties.position.array
    if len(coords)>0 : 
        region = []
        for i in range(3):
            region.append([np.min(coords[:,i])/lat_param,np.max(coords[:,i])/lat_param])
            #print([np.min(coords[:,i])/lat_param,np.max(coords[:,i])/lat_param])
        region = np.array(region)
        if np.max(np.ceil(abs(region))) > (boxL-1):
            print('Warning: cascade crosses simulation box boundary')
            return None #exceed boundary, invalid
        ####np.savetxt("damage_region.txt", np.sign(region)*np.ceil(abs(region)),fmt='%d')
    
    
    # cluster analysis
    # vacancy cluster cutoff (NN4)  sqrt(11)/2*a
    node.modifiers.append(ClusterAnalysisModifier(cutoff = np.sqrt(11.0)/2.0*lat_param, sort_by_size = True))
    # SIA cluster cutoff (NN3) sqrt(2)*a; 
    #node.modifiers.append(ClusterAnalysisModifier(cutoff = np.sqrt(2.0)*lat_param, sort_by_size = True))
    node.compute()
    cluster_sizes = np.bincount(node.output.particle_properties['Cluster'].array)
    ####print('Number of cluster: %d' %(len(cluster_sizes)-1))#excluding cluster ID 0
    ####print('Total defects : %d' %np.sum(cluster_sizes))
    ####np.savetxt("cluster_sizes.txt", cluster_sizes[1:],fmt='%d')
    return cluster_sizes[1:]


if __name__ == '__main__':
    cluster_sizes = []
    N = 20 #N simulations
    reference_filename = 'dump_minimized.data'
    boxL = int(sys.argv[1])
    effective_simulations = 0
    for i in range(N):
        filename = 'dump_irra_minimized_'+str(i)+'.data'
        cluster = cluster_analysis(filename,reference_filename,boxL)
        if cluster is not None: 
            effective_simulations += 1
            cluster_sizes.append(list(cluster))
    cluster_sizes = list(itertools.chain.from_iterable(cluster_sizes)) #combine all results together
    counter = collections.Counter(cluster_sizes) #frequency
    cluster_freq = [] 
    cluster_freq.append(np.array(list(counter.keys())))
    cluster_freq.append(np.array(list(counter.values()))*1.0/effective_simulations) #average frequency
    print('effective simulations: ',  effective_simulations)
    cluster_freq = np.transpose(np.array(cluster_freq))
    cluster_freq = cluster_freq[cluster_freq[:,0].argsort()]
    np.savetxt('vcluster_frequency.txt',cluster_freq,fmt='%d %.4f')
