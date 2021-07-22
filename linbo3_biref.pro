function linbo3_biref,wave,temp

;  procedure to compute the birefringence of linbo3
;  from work of Schlarb and Betzler, Phys Rev B, 48,
;  Num 21, 15613, 1993. The calculation is valid for
;  temperatures of 50-600K and wavelengths of 400-1200 nm

;  the temperature (C) and wavelength (nm) are input

n_e=linbo3_eindex(wave,temp)
n_o=linbo3_oindex(wave,temp)

return,n_e-n_o

end