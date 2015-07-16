''' 
Program to calculate Cosmic Variance in 6dFGSv Zeropoint
'''

import numpy as np
import cosmologycodes as cc
from astrofunctions import getRA, getDec
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

def MLE_BF(g):
    ngal = len(g)
    x_arr = g[:,0]
    y_arr = g[:,1]
    z_arr = g[:,2]
    r_arr = g[:,3]
    v_arr = g[:,6]
    sigmastar2 = 250.**2
    sigma_n = np.zeros(ngal)

    Aij = np.zeros([3,3])
    rihat = np.array([x_arr/r_arr, y_arr/r_arr, z_arr/r_arr])
    for i in range(3):
        for j in range(3):
            Aij[i,j] = np.sum(rihat[i,n]*rihat[j,n] / (sigma_n[n]**2. + sigmastar2) for n in range(ngal))
    Aij_inv = np.linalg.pinv(Aij)
    rsum = np.array([np.sum(rihat[0,n]*v_arr[n]/(sigma_n[n]**2.+sigmastar2) for n in range(ngal)),
                     np.sum(rihat[1,n]*v_arr[n]/(sigma_n[n]**2.+sigmastar2) for n in range(ngal)),
                     np.sum(rihat[2,n]*v_arr[n]/(sigma_n[n]**2.+sigmastar2) for n in range(ngal))])
    u = np.dot(Aij_inv,rsum)
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

def meanVector(tuplelist):
    meanvec = [np.mean(y) for y in zip(*tuplelist)]
    return meanvec

def linearVelocityTest(mocksFolder):
    ''' Test Normalising Velocities to Zero '''

    vArr = []
    magArr = []
    BFdiffArr = []
    dmagArr = []
    dDecArr = []
    dRAArr = []
    # Read in each mock. 
    for i in range(20):
        mockFile = '%svmock_Gz_highbias_c0p6_%s.dat' % (mocksFolder,str(i).zfill(2))
        data = np.loadtxt(mockFile, comments='#')
        g = data[:,[3,4,5,6,7,8,11]]
        # Calculate true SSBF
        trueBF = MLE_BF(g)
        trueDec = getDec(trueBF[0],trueBF[1],trueBF[2])
        trueRA = getRA(trueBF[0],trueBF[1])

        # Find galaxies within -20 < Dec < 0
        gInC = g[g[:,5]>-20]
        # Calculate mean radial velocity within -20 < Dec < 0 circle
        v = np.mean(gInC[:,6])
        # Renormalise all velocities in mock to normalise Great Circle to zero
        g[:,6] -= v
        # New SSBF
        newBF = MLE_BF(g)
        newDec = getDec(newBF[0],newBF[1],newBF[2])
        newRA = getRA(newBF[0],newBF[1])

        BFdiff = newBF-trueBF
        dDec = newDec - trueDec
        dRA = newRA - trueRA

        #print 'Old BF: ', trueBF, ', New BF: ', newBF
        print 'Change in BF: ', mag(BFdiff)
        print 'Old Dec: ', trueDec, 'New Dec: ', newDec
        print 'Change in Dec: ', dDec

        dmag = mag(newBF) - mag(trueBF)
        dmagArr.append(dmag)
        dDecArr.append(dDec)
        dRAArr.append(dRA)

        magArr.append(mag(BFdiff))
        BFdiffArr.append((BFdiff[0],BFdiff[1],BFdiff[2]))
        vArr.append(v)

    vArr = np.array(vArr)
    magArr = np.array(magArr)
    dmagArr = np.array(dmagArr)
    dDecArr = np.array(dDecArr)
    dRAArr = np.array(dRAArr)

    meanBFdiff = meanVector(BFdiffArr)

    print ''
    print 'Mean radial v in circle: ', np.mean(vArr)
    print 'Variance: ',np.var(vArr)
    print 'Std: ', np.std(vArr)
    print ''

    print 'Mean BF diff mag: ', np.mean(magArr)
    print 'Std: ', np.std(magArr)
    print ''

    print 'Mean BF diff: ', meanBFdiff
    print 'Mag: ', mag(meanBFdiff)
    print ''

    print 'Mean dmag: ', np.mean(dmagArr)
    print 'Std: ', np.std(dmagArr)

    print 'Mean dDec: ', np.mean(dDecArr)
    print 'Std: ', np.std(dDecArr)

    print 'Mean dRA: ', np.mean(dRAArr)
    print 'Std: ', np.std(dRAArr)

    return

def observableEtaTest(mocksFolder):
    ''' Test Normalising Observed eta to Zero '''

    zfile = '/Users/Morag/Cosmology_codes/redshift_OM0.3175.dat'

    xArr = []
    magArr = []
    BFdiffArr = []
    # Read in each mock. 
    #  r  s  i  x[Mpc/h]  y[Mpc/h]  z[Mpc/h]  D_c[Mpc/h]  RA  Dec  ztrue  zobs  vpec[km/s]  appJ  V_max_sub
    for i in range(20):
        mockFile = '%svmock_Gz_highbias_c0p6_%s.dat' % (mocksFolder,str(i).zfill(2))
        data = np.loadtxt(mockFile, comments='#')
        # x[Mpc/h]  y[Mpc/h]  z[Mpc/h]  D_c[Mpc/h]  RA  Dec vpec[km/s] zobs
        g = data[:,[3,4,5,6,7,8,11,10]]
        # Calculate true SSBF
        trueBF = MLE_BF(g)

        # Calulate observable x 
        x = v_to_x(g[:,6], g[:,7], g[:,3])
        # Append x values to g
        g = np.column_stack((g, x[:,None]))

        # Find galaxies within -20 < Dec < 0
        gInC = g[g[:,5]>-20]
        # Calculate mean radial velocity within -20 < Dec < 0 circle
        meanx = np.mean(gInC[:,8])
        
        # Renormalise all x values in mock to normalise Great Circle to zero
        g[:,8] -= meanx

        # Convert these back to velocities
        OM = 0.31
        Dz = cc.dxx_mpch_arr(len(g), g[:,7], OM)
        g[:,6] = x_to_v(g[:,8], g[:,7], Dz, zfile)

        # Calculate true SSBF
        newBF = MLE_BF(g)

        BFdiff = newBF-trueBF
        print 'Old BF: ', trueBF, ', New BF: ', newBF
        print 'Change in BF: ', mag(BFdiff)
        magArr.append(mag(BFdiff))
        BFdiffArr.append((BFdiff[0],BFdiff[1],BFdiff[2]))
        xArr.append(meanx)

    xArr = np.array(xArr)
    magArr = np.array(magArr)

    meanBFdiff = meanVector(BFdiffArr)

    print 'Mean radial x in circle: ', np.mean(xArr)
    print 'Variance: ',np.var(xArr)
    print 'Std: ', np.std(xArr)
    print ''

    print 'Mean BF diff mag: ', np.mean(magArr)
    print 'Std: ', np.std(magArr)
    print ''

    print 'Mean BF diff: ', meanBFdiff
    print 'Mag: ', mag(meanBFdiff)
        
    return


def main():

    mocksFolder = '/Users/Morag/6dFGS/Mock_code/mocks/'

    ''' Test Normalising Velocities to Zero '''
    linearVelocityTest(mocksFolder)

    ''' Test Normalising Observed eta to Zero '''
    #observableEtaTest(mocksFolder)



if __name__ == "__main__":
    main()