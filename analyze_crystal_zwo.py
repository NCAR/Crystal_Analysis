import os
import numpy
import glob

global fit, xfit, yfit

#Setup
debug = "no"
ans = ' '
in_dir = 'E:\\Crystal Test\\'
out_dir = 'C:\\Users\\tomczyk\\Documents\\HAO-IG\\Crystal Evaluation\\Crystal Scan\\Results\\'

os.chdir('C:\Users\tomczyk\Documents\HAO-IG\Crystal Evaluation\Crystal Scan')

#options TODO Maybe break these out into an ini file
rem_poly=False  #remove polinomial Fit True/False
np = 4 # number of polinomials to be removed
rem_zernike = False     #remove zernike True/False
nz = 35 # remove zernike polinomial upto 35
rem_trend=True # remove first order trend
if rem_trend:
	np = 1 # if the above is set it should set this to 1

#TODO original code is 632.d is this notation just to specify that it should be a float rather than an int
wave0 = 632.0
npixel = 4640				   #number of CCD columns, ZWO camera

#TODO original code reads clear = 0.8d
clear = 0.8          #clear aperture fraction

#TODO why is nfiles specified as a double wont it always be an int?
files = glob.glob(in_dir+'202006*.txt')
nfiles = double(len(files))

#Iterate through Files and open

for file in files:
    print(file)
    infile = open(file,"r")
    desc=infile[0]
    #TODO need to find what to split at
    thickness,in_crystal = infile[1].split()
    nx, ny, nspec, diam = infile[2].split()

    #format crystal string
    crystal = in_crystal.replace(" ", "").lower()

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

    bin = npixel/nspec #CCD binning
    print("bin:", bin)

# select dispersion for grating used
    dispersion_dict = {
        '2400': 0.00185*bin,  #dispersion, (nm/pixel) 2400 l/mm grating
        '1800': 0.00214*bin,
        '600': 0.01144*bin,
        'echelle': 0.00119*bin
    }
    disp=dispersion_dict.get(grating)
    #TODO replaced IDL dindgen with numpy arrange need to test to make sure its working as expected
    wave = numpy.arrange(nspec) * disp + wave0 - nspec * disp / 2
    #TODO need to convert the calcite brief and other functions to python
	#TODO Could not find quartz_biref file is this not used?
    #case crystal of
	# 	'calcite': begin
	# 					mu=calcite_biref(wave0,25.d)
	# 					dlam = 10.d
	# 					dmu = calcite_biref(wave0+dlam/2.d,25.d) - calcite_biref(wave0-dlam/2.d,25.d)
	# 					dmu_dlam = dmu/dlam         ;derivative of birefringence (nm-1)
	# 					eff_mu=mu*(dmu_dlam*wave0/mu -1.d)			;effective birefringence (see Title, SoPh, 39, 505.)
	# 					fsr=wave0^2/abs(eff_mu*thickness*1.d6)
	# 					ld0=eff_mu*thickness*1.d6
	# 				end
    #
	# 	'linb':	begin
	# 					mu=linbo3_biref(wave0,25.d)
	# 					mu2=linbo3_biref(wave0+1.d,25.d)
	# 					dmu_dlam=mu2-mu     ;take derivative of birefringence (nm-1)
	# 					eff_mu=mu*(dmu_dlam*wave0/mu -1.d)			;effective birefringence (see Title, SoPh, 39, 505.)
	# 					fsr=wave0^2/abs(eff_mu*thickness*1.d6)
	# 					ld0=eff_mu*thickness*1.d6
	# 				end
    #
	# 	'quartz': begin
	# 					mu=quartz_biref(wave0,25.d)
	# 					mu2=quartz_biref(wave0+1.d,25.d)
	# 					dmu_dlam=mu2-mu     ;take derivative of birefringence (nm-1)
	# 					eff_mu=mu*(dmu_dlam*wave0/mu -1.d)			;effective birefringence (see Title, SoPh, 39, 505.)
	# 					fsr=wave0^2/abs(eff_mu*thickness*1.d6)
	# 					ld0=eff_mu*thickness*1.d6
	# 				end
    #
	# endcase

    ifsr = int(fsr/disp)		#convet free spectral range to pixel
    if debug:
        print(fsr, ifsr)
    #TODO Confirm that dblarr is analogous to numpy.zeros
    ldelta = numpy.zeros(nx, ny)
    intens = numpy.zeros(nx, ny)
    x = numpy.zeros(nx, ny)
    y = numpy.zeros(nx, ny)


    # read Transmission Spectra
#TODO Find analogous to window maybe just set up for a graphic matplotlib
#  window,0,xs=690,ys=500,xpos=0,ypos=0
#  window,6,xs=690,ys=500,xpos=0,ypos=0

    trans=numpy.zeros(nspec)
    #TODO 3 Dimensional Array? data=dblarr(nx,ny,nspec,/nozero)
    data=numpy.zeros(nx, ny, nspec)

#
# ;  file format:
# ;
# ;  line 1: descriptor
# ;  line 2: crystal thickness (mm), crystal type
# ;  line 3: x points in scan, y points in scan, number of points in spectrum, crystal diameter (mm)
# ;  then there are nx by ny lines with the first 2 numbers being the x position and y position followed
# ;  by the spectrum

def fit_lyot(p):
    model = p[0]*numpy.cos(numpy.pi*p[1]/xfit)**2
    chisq=sum((model-yfit)**2)
    return chisq

#Returns an array?
def ldfunc(params):
    fit, xfit, yfit
    x = (xfit-min(xfit))/(max(xfit) - min(xfit))
    f = params[0] + params[1] * x + params[2] / (x - params[3]) + params[4] * x ^ 3
    return f

#returns a number?
def ldmin(params):
    global fit, xfit, yfit
    return sum((ldfunc(params)-yfit)**2)
