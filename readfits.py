import math
import pyfits
import numpy as np
import glob,os
import matplotlib.pyplot as plt


#W3IRS5
#1-flux1 //2116
#2-wave1 //2116
#3-flux2 //2116
#4-flux3 //2116
#6-ra    //5x5
#7-ra   //5x5
#8-dec   //5x5
#9-dec   //5x5
#11-flux4  //2116
#12-wave   //2116 -> exactly the same as 2
#14-flux5  //1058 -> looks like 2-binned channel 1, but it's not
#15-wave5  //1058 -> looks like 2-binned channel 2, but it's not
#16-flux6  //1058 -> 3
#17-flux7  //1058 -> 4
#19-ra   //5x5 -> perfectly the same as 6
#20-ra   //5x5
#21-dec  //5x5
#22-dec  //5x5
#24-flux //1058 -> 11 
#25-wave //1058 -> 2

#ok, we have two sets of data, one having double spectral
#resolution
#flux is of course flattened at low resolution, which may result
#in flux differences

#each set contains four sets of fluxes with similar shape
#but a bit different values




def fluxtest(flux,wave):
	fluxtab = np.empty(len(flux))
	wavetab = np.empty(len(flux))
	lol = 0
	for j in range(len(wave)):
		if math.isnan(flux[j,0,0])==False:
			fluxtab[lol] = flux[j,0,0]
			wavetab[lol] = wave[j]
			lol = lol+1
	fluxtab=fluxtab[0:lol]
	wavetab=wavetab[0:lol]
	plt.plot(wavetab,flu1xtab)
	return





fitslist=glob.glob('*.fits') #search current catalog for fits files
print fitslist
for i in fitslist:
	hdulist = pyfits.open(i)	

	header = hdulist[0].header	#header
	ra1 = hdulist[6].data
	ra2 = hdulist[7].data
	dec1 = hdulist[8].data
	dec2 = hdulist[9].data
		
#	print 'flux1', len(hdulist[1].data)
#	print 'wavelength', len(hdulist[2].data)
#	print 'flux2', len(hdulist[3].data)
#	print 'flux3', len(hdulist[4].data)
#	print '?', len(hdulist[5].data)
#	print 'ra1', len(hdulist[6].data)
#	print 'dec1', len(hdulist[7].data)
#	print 'ra2', len(hdulist[8].data)
#	prin1t 'dec2', len(hdulist[9].data)


#	print header
	print hdulist.info()
	flux = hdulist[1].data[:,2,2]
	wave = hdulist[2].data
#	not_nan = np.where(np.isnan(flux) == False)
#	flux = flux[not_nan]
#	wave = wave[not_nan]
#	lines = hdulist[5].data[not_nan]


#	ind0 = np.where(lines == 9)
#	flux0 = flux[ind0]
#	wave0 = wave[ind0]
#
#	print flux0, wave0

#	print np.isnan(hdulist[1].data)
#	index = np.where(hdulist[5].data == 9)[0]
#	print hdulist[5].data[index]
#	print hdulist[2].data[index]
#	print hdulist[2].data
#	print hdulist[1].data[index,0,0]

	plt.plot(wave,flux,'o')


#	fluxtest(hdulist[3].data,hdulist[2].data)
#	fluxtest(hdulist[4].data,hdulist[2].data)
#	fluxtest(hdulist[11].data,hdulist[12].data)
#	fluxtest(hdulist[14].data,hdulist[15].data)
#	fluxtest(hdulist[16].data,hdulist[15].data)
#	fluxtest(hdulist[17].data,hdulist[15].data)
#	fluxtest(hdulist[24].data,hdulist[25].data)
	
	plt.savefig(i+'test.eps',format='eps',dpi=100)
	plt.clf()
