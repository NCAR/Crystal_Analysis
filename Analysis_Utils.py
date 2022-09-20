#     file format:
# ;  
# ;  line 1: descriptor
# ;  line 2: crystal thickness (mm), crystal type
# ;  line 3: x points in scan, y points in scan, number of points in spectrum, crystal diameter (mm)
# ;  then there are nx by ny lines with the first 2 numbers being the x position and y position followed
# ;  by the spectrum 
# read the first line in the infile and put it into descriptor
def read_header(filename):
    with open(filename, 'r') as f:
        desc = f.readline()
        thickness, in_crystal = f.readline().split()
        nx, ny, nspec, diam = f.readline().split()
        crystal = in_crystal.replace(" ", "").lower()
        return desc, thickness, crystal, nx, ny, nspec, diam

#Get the Dispersion of a grating
def get_dispersion(grating):
    if grating == '2400':
        disp = 0.00185*bin
    elif grating == '1800':
        disp = 0.00214*bin
    elif grating == '600':
        disp = 0.01144*bin
    elif grating == 'echelle':
        disp = 0.00119*bin
    return disp