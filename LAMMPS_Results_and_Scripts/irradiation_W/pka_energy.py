import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=24)


class MDPKAEnergy:
    def __init__(self):
        self.A = 183.84
        self.Z = 74

    #calculate the energy used to displace atoms
    #according to Was book P109/839 equation (2.54)
    #natural units
    #energy_arr should be in eV
    def MDEnergy(self,energy_arr):
        g = lambda epsilonN: 3.4008* epsilonN**(1.0/6)+0.40244*epsilonN**(3.0/4)+epsilonN
        k = lambda Z1,Z2,A1: 0.1337*Z1**(1.0/6)*(Z1/A1)**0.5
        #T in eV
        epsilon = lambda Z1,Z2,A1,A2,a,T: (A2*T/(A1+A2))*a/(Z1*Z2*8.5424546e-2**2)
        #bohr radius in 1/eV
        radius = lambda Z1,Z2: (9.0*np.pi**2/128)**(1.0/3)*2.6817268e-4*(Z1**(2.0/3)+Z2**(2.0/3))**(-0.5)
        A1 = A2 = self.A
        Z1 = Z2 = self.Z
        Emd = np.zeros(energy_arr.shape)
        for i,T in enumerate(energy_arr):
            a = radius(Z1,Z2)
            epsilonN = epsilon(Z1,Z2,A1,A2,a,T)
            Emd[i] = T/(1.0+k(Z1,Z2,A1)*g(epsilonN))
        return Emd

if __name__ == '__main__':
    
    ##ion irradiation
    #energies = [0.1,0.16,0.25,0.4,0.63,1, 1.58,2.51,4, 6.3,10,15.8, 25.1,40,63, 100, 150, 250, 400] #in keV
    data = np.genfromtxt('pka_150keV.txt',skip_header=1)
    energies = list(data[:,2])
    #energies = np.loadtxt('tmp.txt').flatten()


    ##neutron irradiation
    #energies = [0.15,0.25,0.4,0.75,1,1.5,2.51,4,5,7.5,10,12.6,15.8,25.1,30,40,60,75,100,150]
    #energies = np.loadtxt('tmp.txt').flatten()

    PKA = MDPKAEnergy()
    energies = np.array(energies)*1.0e3 #in eV
    Emd = PKA.MDEnergy(energies)

    print 'Epka vs MD energy in keV: '
    for i in range(len(energies)):
        print Emd[i]/1000.0
