import numpy as np
import sys


#command: e.g. python generate_random_directions.py 20 0.1
def generateRandomDirections(N,halfL,Energy):
    m = 3.0527348e-25 #mass in kg
    E = Energy*1000.0*1.602e-19 #energy in J
    v = np.sqrt(2.0*E/m) /100.0 #velocity in A/ps
    data = np.zeros((N,4))
    data[:,0] = halfL #half length of simulation box
    A = np.array([0.0,1.0,0.0])
    B = np.array([-0.5,0.5,0.0])
    C = np.array([-0.5,0.5,0.5])
    np.random.seed(173257393)
    for i in range(N):
       rand = np.random.rand(2)
       P = (1.0-np.sqrt(rand[0]))*A+np.sqrt(rand[0]*(1.0-rand[1]))*B+rand[1]*np.sqrt(rand[0])*C 
       #ref : http://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
       P = P/np.linalg.norm(P)*v
       data[i,1:] = np.copy(P)
    np.savetxt('velocities.txt',data,fmt='%d %.4f %.4f %.4f')

def generateSpeed(delta=0.5):
    Energy = np.array([30,40,50,60])/1000.0 #keV
    m = 3.0527348e-25 #mass in kg
    E = Energy*1000.0*1.602e-19 #energy in J
    v = np.sqrt(2.0*E/m) /100.0 #velocity in A/ps
    r1 = np.arange(0.5,1.5,0.1)
    r2 = r1 + delta
    data = np.zeros((len(E)*len(r1),4))
    n = 0
    for j in range(len(E)):
        for i in range(len(r1)):
            data[n,0] = Energy[j]
            data[n,1] = v[j]
            data[n,2] = r1[i]
            data[n,3] = r2[i]
            n = n + 1
    np.savetxt('study_ZBL.txt',data,fmt='%.4f %.4f %.4f %.4f')


if __name__ == '__main__':
    #sys.argv[1]: half length of simulation box
    #sys.argv[2]: energy in keV
    if len(sys.argv) != 3 :
        raise Exception('Argument number not correct')
    generateRandomDirections(30,int(sys.argv[1]),float(sys.argv[2]))
    #generateSpeed(delta=0.5)
