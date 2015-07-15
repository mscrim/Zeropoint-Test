import numpy as np

def getDec(x, y, z):
	r = np.sqrt(x**2 + y**2 + z**2)
	Dec = 90. - acos(z/r) * 180. / pi
	return Dec

def getRA(x,y):
	RA = atan2(y,x) * 180./pi
	if RA < 0: RA += 360.
	return RA
