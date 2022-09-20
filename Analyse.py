import os
import glob
from Analysis_Utils import *
import math

global fit, xfit, yfit
#User Defined Variables
debug = False
rem_poly = False
rem_zernike = False
rem_trend = True
in_dir = 'E:\Crystal Test\'
out_dir = 'C:\Users\tomczyk\Documents\HAO-IG\Crystal Evaluation\Crystal Scan\Results\'
#cd into a directory
os.chdir('C:\Users\tomczyk\Documents\HAO-IG\Crystal Evaluation\Crystal Scan')


ans = ''

#Remove Polynomial fit (True or False)
if rem_poly:
    np = 4 #Order Poluyomial to be removed

#Remove Zernike fit (True or False)
if rem_zernike:
    nz = 35 #Remove Zernike polynomial up to 35

#Remove first order trend (True or False)
if rem_trend:
    np = 1 #If the above is set it should set this to 1

#Centeral Wavelength
wave0 = 632.0
#Number of CCD columns, ZWO camera
npixel = 4640
#Clear aperture fraction
clear = 0.8

#Find files to analyze
files = glob.glob(in_dir+'202006*.txt')
nfiles = len(files)

for i in range(nfiles):
    infile = files[i]
    print(infile)

    #read in the Header to the file
    desc, thickness, crystal, nx, ny, nspec, diam = read_header(infile)

    #set grating
    grating = '1800' # Grating used is ('2400', '600', '1800', or 'echelle')
    if '600' in desc:
        if debug:
            print(desc)
            print(thickness)
            print('crystal:'+crystal)
            print(nx, ny, nspec)

    #Setup for plotting
    factor = int(400./nx+0.5)
    factor = 10
    print("factor:"+factor)

    #CCD binning
    bin = npixel/nspec
    print(f"bin: {bin}")

    #Select dispersion fo grating
    disp = get_dispersion(grating)
    #Create an array from 0 to nspec-1 with a step of 1 using double precision
    wave = np.arange(0,nspec,1,dtype=np.float64)*disp+wave0-disp*nspec/2.

    #Nominal free spectral rance 1d nm
    #Caser Switch for Type of Crystal
    if crystal == 'calcite':
        mu = calcite_biref(wave, 25.0)
        dlam = 10.
        dmu = calcite_biref(wave0+dlam/2.0,25.0)-calcite_biref(wave0-dlam/2.0,25.0) 
        dmu_dlam = dmu/dlam         #derivative of birefringence (nm-1)
        eff_mu= mu * (dmu_dlam*wave0/mu-1.0)	#effective birefringence (see Title, SoPh, 39, 505.)
        fsr= wave0^2 / abs(eff_mu*thickness*1.0)
        ld0=eff_mu*thickness*1.0
    
    if crystal == 'linb':
        mu = linbo3_biref(wave, 25.0)
        mu2 = linbo3_biref(wave0+1.0,25.0)
        dmu = mu2-mu #derivative of birefringence (nm-1)
        eff_mu = mu*(dmu_dlam*wave0/mu-1.0)    #effective birefringence (see Title, SoPh, 39, 505.)
        fsr= wave0^2 / abs(eff_mu*thickness*1.0)
        ld0=eff_mu*thickness*1.0
    
    if crystal == 'quartz':
        mu = quartz_biref(wave, 25.0)
        mu2 = quartz_biref(wave0+1.0,25.0)
        dmu = mu2-mu #derivative of birefringence (nm-1)
        eff_mu = mu*(dmu_dlam*wave0/mu-1.0)    #effective birefringence (see Title, SoPh, 39, 505.)
        fsr= wave0^2 / abs(eff_mu*thickness*1.0)
        ld0=eff_mu*thickness*1.0

    ifsr = fix(fsr/disp) #Convert free spectral range to pixels

    if debug: print(f"{fsr},{ifsr}")
    #Create a nx by ny zero array for the data
    ldelta = np.zeros((nx,ny),dtype=np.float64)
    intens = np.zeros((nx,ny),dtype=np.float64)
    x=np.zeros((nx,ny),dtype=np.float64)
    y=np.zeros((nx,ny),dtype=np.float64)

    #Read in transmission spectra
    window,0,xs=690,ys=500,xpos=0,ypos=0
    window,6,xs=690,ys=500,xpos=0,ypos=0

    trans = np.zeros((nspec),dtype=np.float64)
    data = np.zeros((nx, ny, nspec))
    
    #Read in all the data from the infile line by line and place them in the data array
    for iy in range(ny):
        for ix in range(nx):
            #Read the next line in infile and assign the three comma seperated valuse to xpos and ypos and trans
            xpos, ypos, trans = infile.readline().split()

            data[ix,ny-1-iy,:]=trans
            intens[ix,ny-1-iy]=max(trans)-min(trans)
            x[ix,iy]=xpos
            y[ix,iy]=ypos

            if debug: print(f"{ix},{iy},{xpos},{ypos}")
    dx = abs(x[1,0]-x[0,0])
    dy = abs(y[1,0]-y[0,0])

    x=x-x[0,0]
    y=y-y[0,0]

    #Fit the Transmission Spectra
    #Obtain guess for ldelta using grid search applied to spectrum near center of crystal
    num = 10000
    chisq = np.zeros((num),dtype=np.float64)
    int_thresh = 0.5 #intensity Threshold
    ld = ld0*(0.85+0.3*(np.arange(0,num,1,dtype=np.float64)/double(num-1)))

    trans = data[nx/2, ny/2, :] #Get Central Spectrum
    
    #TODO Fix this
    wset,6
    plot, wave, trans

    rmin = run_min(trans,ifsr*1.0)
    rmax = run_max(trans,ifsr*1.0)
    trans = (trans-rmin)/(rmax-rmin)

    for j in range(num-1):
	    model = cos(!pi*ld[j]/wave)^2
		chisq[j] = total( (model-trans)^2 )

    mn = min(cisq, imin)
    guess = [1.0, ld[imin]]

    model = cos(!pi*ld[imin]/wave)^2
    wset,0
    plot,wave,trans,xsty=1,yr=[0,1],ysty=1,psym=3
	oplot,wave,model
	if debug: ans = input.read("enter return")


