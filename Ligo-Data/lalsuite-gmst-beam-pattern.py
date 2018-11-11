import lal

#Python code to convert gpstime to GMST and to find beam pattern functions.
# Install lalsuite by typing "pip install lalsuite" in command line
#To run the script just do 'python lal-gmst-beam-pattern.py'

timeToUse = 1126259462.41302 #GPStime
IFO_l = lal.CachedDetectors[lal.LALDetectorIndexLLODIFF] #Livingston
IFO_h = lal.CachedDetectors[lal.LALDetectorIndexLHODIFF] #Handoford
fancy_time = lal.lal.LIGOTimeGPS(timeToUse) #Nothing important
#print fancy_time
gmst = lal.GreenwichMeanSiderealTime(fancy_time) #converts gpstime to gmst
print gmst
RA = 1.123 #right-ascension
dec = -0.012 #declination
#change the RA and dec to the numbers from the posteior files (the ones I use here are just random...)
psi = 1.57 #polarization
F_plus_l, F_cross_l= lal.ComputeDetAMResponse(IFO_l.response, RA, dec, psi, gmst) #beam pattern for Livingston
F_plus_h, F_cross_h= lal.ComputeDetAMResponse(IFO_h.response, RA, dec, psi, gmst) #beam pattern for Hanford
print F_plus_h, F_cross_h
print F_plus_l, F_cross_l
