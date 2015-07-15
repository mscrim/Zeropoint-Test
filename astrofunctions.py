import numpy as np

def getDec(x, y, z):
	r = np.sqrt(x**2 + y**2 + z**2)
	Dec = 90. - np.arccos(z/r) * 180. / np.pi
	return Dec

def getRA(x,y):
	RA = np.arctan2(y,x) * 180./np.pi
	if RA < 0: RA += 360.
	return RA
