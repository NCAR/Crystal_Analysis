function calcite_biref,wave,t

;  routine to return birefringence of calcite from Beckers and Dunn paper
;  wave is wavelength in nm, t is temperature in degrees C

w=wave/1.d3
mu=-0.163724d0 -3.15d-3/w^2 -3.896d-5/w^4 -2.911d-6/w^6 +3.037d-3*w^2 $
 +2.54d-4*w^4 -2.52d-5*w^6 +1.d-5*(t*(1.044-0.16*w)+0.00043*t^2)

return,mu
end
