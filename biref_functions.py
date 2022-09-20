
# ;  routine to return birefringence of calcite from Beckers and Dunn paper
# ;  wave is wavelength in nm, t is temperature in degrees C
def calcite_biref(wave, t, debug = False):
    w=wave/1.00
    mu=-0.1637240 -3.15-3/w^2 -3.896-5/w^4 -2.911-6/w^6 +3.037-3*w^2+2.54-4*w^4 -2.52-5*w^6 +1.0-5*(t*(1.044-0.16*w)+0.00043*t^2)
    if debug: print(f"calcite_biref: {mu}")
    return mu



# ;  procedure to compute the birefringence of linbo3
# ;  from work of Schlarb and Betzler, Phys Rev B, 48,
# ;  Num 21, 15613, 1993. The calculation is valid for
# ;  temperatures of 50-600K and wavelengths of 400-1200 nm
# ;  the temperature (C) and wavelength (nm) are input
def linbo3_biref(wave, t, debug = False):
    n_e=linbo3_eindex(wave,t)
    n_o=linbo3_oindex(wave,t)
    return n_e-n_o

