import numpy as np
import sys


#command: e.g. python generate_random_directions.py 20 0.1
def generateRandomDirections(N,halfL,Energy):
    m = 3.0527348e-25 #mass in g
    E = Energy*1000.0*1.602e-19 #energy in J
    v = np.sqrt(2.0*E/m) /100.0 #velocity in A/ps
    data = np.zeros((N,4))
    data[:,0] = halfL #half length of simulation box
    A = np.array([0.0,1.0,0.0])
    B = np.array([-0.5,0.5,0.0])
    C = np.array([-0.5,0.5,0.5])
    for i in range(N):
       rand = np.random.rand(2)
       P = (1.0-np.sqrt(rand[0]))*A+np.sqrt(rand[0]*(1.0-rand[1]))*B+rand[1]*np.sqrt(rand[0])*C 
       #ref : http://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
       P = P/np.linalg.norm(P)*v
       data[i,1:] = np.copy(P)
    np.savetxt('velocities.txt',data,fmt='%d %.4f %.4f %.4f')


if __name__ == '__main__':
    #sys.argv[1]: half length of simulation box
    #sys.argv[2]: energy in keV
    if len(sys.argv) != 3 :
        raise Exception('Argument number not correct')
    generateRandomDirections(20,int(sys.argv[1]),float(sys.argv[2]))
