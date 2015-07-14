''' 
Program to calculate Cosmic Variance in 6dFGSv Zeropoint
'''

import numpy as np

def galaxiesInCircle(mockFile):
	''' 
	Read in mock catalogue, selecting columns
	x[Mpc/h]  y[Mpc/h]  z[Mpc/h]  D_c[Mpc/h]  RA  Dec  vpec[km/s]
	and select only galaxies with Dec > -20
	'''
	d = np.loadtxt(mockFile, comments='#')
	x = d[:,[3,4,5,6,7,8,11]]
	gInC = x[x[:,5]>-20]
	return gInC

def circleVelocity(gInC):
	''' Calculate mean radial velocity of galaxies in Great Circle '''
	return np.mean(gInC[:,6])

def main():

	mocksFolder = '/Users/Morag/6dFGS/Mock_code/mocks/'

	vArr = []
	# Read in each mock. Find galaxies within -20 < Dec < 0
	for i in range(20):
		mockFile = '%svmock_Gz_highbias_c0p6_%s.dat' % (mocksFolder,str(i).zfill(2))
		gInC = galaxiesInCircle(mockFile)
		vArr.append(circleVelocity(gInC))
	vArr = np.array(vArr)
	print 'Mean radial v in circle: ', np.mean(vArr)
	print 'Variance: ',np.var(vArr)
	print 'Std: ', np.std(vArr)

if __name__ == "__main__":
    main()