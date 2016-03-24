import math
import pyfits
import numpy as np
import mpfitexpr_simple as mp
import glob,os
import matplotlib.pyplot as plt
import readcol


#=====================================================#
# gaussian function of given parameters               #
def gaussian(x,level,h, peak, sig):
		func = level+h*np.exp(-(x-peak)**2/(2*sig**2))
		return func


c = 299792458 #m/s

output = open('intensdata1.dat','w')
errout = open('errdata1.dat','w')
crdout = open('pacscoord1.dat','w')



compint = np.empty(25)
compcont = np.empty(25)
compinterr = np.empty(25)
compconterr = np.empty(25)
compra = np.empty(25)
compdec = np.empty(25)



#list of lines to measure
linestab = readcol.readcol('../lines.dat')

fitslist=glob.glob('*.fits') #search current catalog for fits files




#loop over files
for f in fitslist:

    #loop over molecules
	for m in range(len(linestab)):
		lab_wvl = float(linestab[m][1]) #central wavelenght
		name = linestab[m][0] #name of the line

		fig = plt.figure() #starting a plot



		hdulist = pyfits.open(f) #open a fits file
		header = hdulist[0].header #load a header



        #loop over spaxels
		spx = 1 #spaxel number count, starting from one to avoid confusion
		
		printcheck = 0

		for i in range(0,5):
			for j in range(0,5):
				ax = fig.add_subplot(5,5,spx)
				#remove NaN and choose spaxels
				index = np.where(np.isnan(hdulist[1].data[:,i,j])==False)
				flux = hdulist[1].data[:,i,j][index]
				wave = hdulist[2].data[index]
				ra1 = hdulist[6].data[:,i,j][index]
				ra2 = hdulist[7].data[:,i,j][index]
				dec1 = hdulist[8].data[:,i,j][index]
				dec2 = hdulist[9].data[:,i,j][index]

				#remove repeating wavelengths
				uniq = np.unique(wave,return_index=True)
				flux = flux[uniq[1]]
				wave = wave[uniq[1]]
			
	
				#subplot
				ax = fig.add_subplot(5,5,spx)


			    #check if this file is of interest for this line
				if (lab_wvl < max(wave) and lab_wvl > min(wave)):
					printcheck = 1

				
					#finding channel closest to lab_wvl
					mintab = np.empty(2000)
					for mm in range(len(wave)):
						mintab[mm] = np.fabs(wave[mm]-lab_wvl)
					mintab = mintab[0:mm]


					lab_index = np.where(mintab == min(mintab))[0][0]

					#plotting only limited area around central line

					cut = 50
					x = wave[lab_index-cut:lab_index+cut]
					y = flux[lab_index-cut:lab_index+cut]

					#channel size



			
					ra1 = ra1[lab_index-cut:lab_index+cut]
					ra2 = ra2[lab_index-cut:lab_index+cut]
					dec1 = dec1[lab_index-cut:lab_index+cut]
					dec2 = dec2[lab_index-cut:lab_index+cut]
					

					#cutting all far channels 
					#for lines on the edge of observing channel

					lim1 = np.where( x > lab_wvl - 0.6)
					lim1_0 = np.where( x[lim1] < lab_wvl + 0.6)

					x = x[lim1][lim1_0]
					y = y[lim1][lim1_0]
					ra1 = ra1[lim1][lim1_0]
					ra2 = ra2[lim1][lim1_0]
					dec1 = dec1[lim1][lim1_0]
					dec2 = dec2[lim1][lim1_0]

					#ax.plot(x,y,drawstyle='steps-mid',color='b',linewidth=0.2)
					

					
					dist = np.empty(200)
					for ii in range(len(x)-1):
						dist[ii] = x[ii+1]-x[ii]
						print dist[ii]
					dist = dist[0:ii]
					chan_size = np.median(dist)
				



					print 'Channel size:', chan_size
					vel_res = c*1e-3*chan_size/lab_wvl		
					print 'Velocity resolution:',vel_res
					vel_limit = 150*lab_wvl/(c*1e-3)
					print '200 km/s limit:', vel_limit
						



					# there is a a chance we entered but there was no line
					#if positive escape this one

					try:
						if not x: 
							printcheck = 0
							continue
					except ValueError:
						twojastara = None

					#limits for rms and linear fit

					lim2 = np.where((x < lab_wvl-vel_limit))
					lim3 = np.where((x > lab_wvl+vel_limit))

					x1 = np.append(x[lim2],x[lim3])
					y1 = np.append(y[lim2],y[lim3])
					
					#same case again, if empty leave
					try:
						if not x1: 
							printcheck = 0
							continue
					except ValueError:
						twojastara = None

					#first linear fit based only on data further away from the line (0.05)
					p = np.poly1d(np.polyfit(x1,y1,1))
					#ax.plot(x1,p(x1),linewidth=0.2)
	
					#substracting linear fit
					y2 = y-p(x)							#y2 - all values, continuum subs
					cont1 = p(lab_wvl)

		#			ax.plot(x,y2,drawstyle='steps-mid',color='red',linewidth=0.2)

					
					#cutting channels away from the line again, now for continuum substracted data

					#x1 is corresponding again
#!!!!!!!!!!!!!!!!!!!
					y3 = np.append(y2[lim2],y2[lim3]) #y3 - values for rms, continuum subs
					x3 = np.append(x[lim2],x[lim3])

					std = np.std(y3)
					print 'std y3:',std
					lim_rms = np.where( np.abs(y3) < std*3 )
					std = np.std(y3[lim_rms])
					print 'std y3[lim]:',std

					p2 = np.poly1d(np.polyfit(x3[lim_rms],y3[lim_rms],1))
					y11 = y2-p2(x)							#y2 - all values, continuum subs
#					ax.plot(x,p2(x),'-b',linewidth=0.2)
					ax.plot(x,y11,drawstyle='steps-mid',color='red',linewidth=0.2)

					y12 = np.append(y11[lim2],y11[lim3])
					y12 = y12[lim_rms]
					std=np.std(y12)


					lim6 = np.where(x>lab_wvl-vel_limit)
					lim7 = np.where(x[lim6]<lab_wvl+vel_limit)
		
			
					y5 = y11[lim6][lim7]
					x5 = x[lim6][lim7]



					#print y3

					
					print 'std y12:', std

					ax.plot([lab_wvl,lab_wvl],[-5,5],'--k',linewidth=0.1)

#channel sizec

				

#rozmiar pola do sumowania (ilosc kanalow)
					len1 = len(y5)
					print 'sum channels:',len1	
					int_std = std		
					int_std = int_std*chan_size*1e-6 #Jy*m
					int_std = int_std*1e-26*c/(1e-6*lab_wvl)**2 #W/m-2
					int_std = int_std*1e-4 #W/cm-2

					print 'not integrated std', int_std

					int_std = int_std*len1
					print 'integrated std', int_std
					
					print y5
					print chan_size

					sum1 = sum(y5)
					sum1 = sum1*chan_size
					sum1 = sum1*1e-6 # Jy*m
					print 'sum',sum1
					#zamiana na w/m-2
					sum1 = sum1*1e-26*c/(1e-6*lab_wvl)**2 #W/m-2
					sum1 = sum1*1e-4 #W/cm-2

					
					ax.plot((x[0],x[len(x)-1]),[std*3,std*3],'--r',linewidth=0.2)





					if sum1 > int_std*3: ax.plot(x1[lim_rms],y3[lim_rms],'.r',markersize=0.3)
					else: ax.plot(x1[lim_rms],y3[lim_rms],'.k',markersize=0.3)

					ax.plot(x5,y5,'.r',markersize=0.3)

					

					plt.xticks(fontsize=0.1)
					plt.yticks(fontsize=0.1)
	
						

				#	print y3[lim_rms]

#trzeba wyznaczyc rms z poprawionych wartosci, niekoniecznie dopasowywac continuum jeszcze raz????




				#RMS calculation on cont removed rms only data

					print 'STD: ',std
					
				# line with 3*sigma
				#plt.axhline(y=std*3)
				#plt.axhline()

#NIE MOZESZ TAK ROBIC TU MUSZA BYC OBA WARUNKI NAARAZZZZ

		

#we need to have the same radec table as x4,y4 if we want to use the same indices

#					compra = np.empty(25)
#					compdec = np.empty(25)
#					low = lab_index - 3
#					high = lab_index + 3
#
#					for lst in range(low,high):
#						ratab[k] = ra1[lst]
#						dectab[k] = ra2[lst]
#						k += 1
#					ratab = ratab[0:k]
#					dectab = dectab[0:k]
#					compra[spx] = np.average(ratab)
#					compdec[spx] = np.average(dectab)

		

#przy sumowaniu konieczne wymnozenie przez szerokosc kanalu


#blad mnozony tez przez szerokosc kanalu + ilosc kanalow


					#std to blad kontinuum!!!
								
					
					
					compint[spx-1] = sum1
					compcont[spx-1] = cont1
					compinterr[spx-1] = int_std
					compconterr[spx-1] = std

					bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.5)
	
			
					ax.text(0.1,0.95,str(sum1), fontsize=3,transform=ax.transAxes,bbox=bbox_props)
					ax.text(0.1,0.05,str(3*int_std), fontsize=3,transform=ax.transAxes,bbox=bbox_props)

					

					compra[spx-1] = (np.average(ra1)+np.average(ra2))/2.
					compdec[spx-1] = (np.average(dec1)+np.average(dec2))/2.

				
			#		compra[spx-1] = np.average(np.average(ra1),np.average(ra2))
			#		compdec[spx-1] = np.average(np.average(dec1),np.average(dec2))


					
			#		func = 'p[0]+p[1]*numpy.exp(-(x-p[2])**2/(2*p[3]**2))'
			#		err = np.empty(len(x4))
			#		err[:] = 0.1
			#		results,yfit = mp.mpfitexpr_simple(func,x4,y4,err,[0.1,13,lab_wvl,0.1])
			#		print results
			#		ax.set_xlim([lab_wvl-low,lab_wvl+high])
			#		ax.plot(x4,gaussian(x4,results[0],results[1],results[2],results[3]),'k')
					plt.xticks(fontsize=3)
					plt.yticks(fontsize=3)


					
			#		print name, results[0]

					fig.savefig(f[6:16]+name+'.eps',format='eps',dpi=100)
					spx += 1
		fig.clf()
		if printcheck == 1:
			
			print >> output, name,
			print >> output, "%6.3f" %lab_wvl,
			for ff in range(24):
				print >> output,"%10.8e" % compint[ff],
			print >> output,"%10.8e" %compint[24]
	
			print >> output, "cont%s"%name,
			print >> output, "%6.3f" %lab_wvl,
			for ff in range(24):
				print >> output,"%10.8e" % compcont[ff],
			print >> output,"%10.8e" %compcont[24]

			print >> errout, name,
			for ff in range(24):
				print >> errout,"%10.8e" % compinterr[ff],
			print >> errout,"%10.8e" %compinterr[24]
	
			print >> errout, "cont%s"%name,
			for ff in range(24):
				print >> errout,"%10.8e" % compconterr[ff],
			print >> errout,"%10.8e" %compconterr[24]

			print >> crdout, name,
			for ff in range(24):
				print >> crdout,"%10.8f" % compra[ff],
			print >> crdout,"%10.8f" %compra[24],
			for ff in range(24):
				print >> crdout,"%10.8f" % compdec[ff],
			print >> crdout,"%10.8f" %compdec[24]
			
			print >> crdout, 'cont'+name,
			for ff in range(24):
				print >> crdout,"%10.8f" % compra[ff],
			print >> crdout,"%10.8f" %compra[24],
			for ff in range(24):
				print >> crdout,"%10.8f" % compdec[ff],
			print >> crdout,"%10.8f" %compdec[24]			
		#	fig2 = plt.figure()
		#	for ff in range(25):
		#		if ff == 5: plt.plot(compra[ff],compdec[ff],'or',markersize=2)
		#		elif ff == 6: plt.plot(compra[ff],compdec[ff],'ob',markersize=2)
		#		else: plt.plot(compra[ff],compdec[ff],'ok',markersize=2)
		#	fig2.savefig(f[6:16]+name+'coords.eps',format='eps',dpi=100)
		#	fig2.clf()
output.close()
errout.close()
crdout.close()
	
