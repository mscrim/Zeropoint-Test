''' 
Program to calculate Cosmic Variance in 6dFGSv Zeropoint
'''

import numpy as np
import cosmologycodes as cc
speedlight = 299792.458

def galaxiesInCircle(mockFile):
    ''' 
    Read in mock catalogue, selecting columns
    x[Mpc/h]  y[Mpc/h]  z[Mpc/h]  D_c[Mpc/h]  RA  Dec  vpec[km/s]
    and select only galaxies with Dec > -20
    '''
    data = np.loadtxt(mockFile, comments='#')
    g = data[:,[3,4,5,6,7,8,11]]
    gInC = g[g[:,5]>-20]
    return gInC

def SS_BF(g):
    ''' Simple Sum Bulk Flow '''
    ngal = float(len(g))
    wx = (g[:,0]/g[:,3])/ngal
    wy = (g[:,1]/g[:,3])/ngal
    wz = (g[:,2]/g[:,3])/ngal
    u = np.array([np.sum(wx*g[:,6]),
                  np.sum(wy*g[:,6]),
                  np.sum(wz*g[:,6])])
    return u

def mag(vector):
    if len(vector) != 3: raise Exception('mag function only accepts 3D vector')
    result = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return result

def v_to_x(v_true, zobs, Dr):
    ngal = len(v_true)
    OM = 0.31
    #zr = (1. + z_obs)/(1.+v_true/speedlight) - 1.
    #Dr = cc.dxx_mpch_arr(ngal, zr, OM)
    Dz = cc.dxx_mpch_arr(ngal, zobs, OM)
    x = np.log10(Dz/Dr)
    return x

def x_to_v(x, z_obs, Dz, zfile):
    Dr_der = Dz / 10.**x
    zr_der = cc.z_x_arr(Dr_der, zfile)
    v = speedlight * ((z_obs - zr_der)/(1.+zr_der))
    return v

def linearVelocityTest(mocksFolder):
    ''' Test Normalising Velocities to Zero '''

    vArr = []
    diff = []
    # Read in each mock. 
    for i in range(20):
        mockFile = '%svmock_Gz_highbias_c0p6_%s.dat' % (mocksFolder,str(i).zfill(2))
        data = np.loadtxt(mockFile, comments='#')
        g = data[:,[3,4,5,6,7,8,11]]
        # Calculate true SSBF
        trueBF = SS_BF(g)
        # Find galaxies within -20 < Dec < 0
        gInC = g[g[:,5]>-20]
        # Calculate mean radial velocity within -20 < Dec < 0 circle
        v = np.mean(gInC[:,6])
        # Renormalise all velocities in mock to normalise Great Circle to zero
        g[:,6] -= v
        # New SSBF
        newBF = SS_BF(g)
        print 'Old BF: ', trueBF, ', New BF: ', newBF
        print 'Change in BF: ', mag(newBF - trueBF)
        diff.append(mag(newBF - trueBF))
        vArr.append(v)
    vArr = np.array(vArr)
    diff = np.array(diff)

    print 'Mean radial v in circle: ', np.mean(vArr)
    print 'Variance: ',np.var(vArr)
    print 'Std: ', np.std(vArr)

    print 'Mean BF diff: ', np.mean(diff)
    print 'Std: ', np.std(diff)
    return

def observableEtaTest(mocksFolder):
    ''' Test Normalising Observed eta to Zero '''

    zfile = '/Users/Morag/Cosmology_codes/redshift_OM0.3175.dat'

    vArr = []
    diff = []
    # Read in each mock. 
    for i in range(20):
        mockFile = '%svmock_Gz_highbias_c0p6_%s.dat' % (mocksFolder,str(i).zfill(2))
        data = np.loadtxt(mockFile, comments='#')
        # x[Mpc/h]  y[Mpc/h]  z[Mpc/h]  D_c[Mpc/h]  RA  Dec  zobs vpec[km/s]
        g = data[:,[3,4,5,6,7,8,10,11]]
        # Calculate true SSBF
        trueBF = SS_BF(g)

        # Calulate observable x 
        x = v_to_x(g[:,7], g[:,6], g[:,3])
        # Append x values to g
        g = np.column_stack((g, x[:,None]))

        # Find galaxies within -20 < Dec < 0
        gInC = g[g[:,5]>-20]
        # Calculate mean radial velocity within -20 < Dec < 0 circle
        meanx = np.mean(gInC[:,8])

        # Renormalise all x values in mock to normalise Great Circle to zero
        g[:,8] -= meanx

        # Convert these back to velocities
        g[:,7] = 




        raise Exception('debug')

        #gInC = g[g[:,5]>-20]
        
    return


def main():

    mocksFolder = '/Users/Morag/6dFGS/Mock_code/mocks/'

    ''' Test Normalising Velocities to Zero '''
    #linearVelocityTest(mocksFolder)

    ''' Test Normalising Observed eta to Zero '''
    observableEtaTest(mocksFolder)



if __name__ == "__main__":
    main()