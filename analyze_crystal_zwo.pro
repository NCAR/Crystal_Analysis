pro analyze_crystal_zwo


;  procedure to analyze crystal test data


common fit, xfit,yfit

compile_opt strictarr

debug = 'no'				;debug mode ('yes' or 'no')

ans = ' '
in_dir = 'E:\Crystal Test\'
out_dir = 'C:\Users\tomczyk\Documents\HAO-IG\Crystal Evaluation\Crystal Scan\Results\'

CD, 'C:\Users\tomczyk\Documents\HAO-IG\Crystal Evaluation\Crystal Scan'

;  options

rem_poly = 'no'				;remove polynomial fit? ('yes' or 'no')
  if rem_poly eq 'yes' then np = 4    ;order of polynomial to be removed
  
rem_zernike = 'no'			;remove zernike fit? ('yes' or 'no')
  nz = 35               ;remove Zernike polynomials up to 35

rem_trend = 'yes'				;remove first order trend 
  if rem_trend eq 'yes' then np = 1    ;remove 1st order polynomial trend

wave0 = 632.d			     ;central wavelength
npixel = 4640				   ;number of CCD columns, ZWO camera

clear = 0.8d          ;clear aperture fraction

;  find files to analyze

files = file_search(in_dir+'202006*.txt')
nfiles = double(n_elements(files))


for ifile=0,nfiles-1 do begin

;  open file

	print,files[ifile]
	openr,1,files[ifile]

	desc = ' '
	thickness = 0.d
	diam = 0.d
	in_crystal = ' '

;  file format:
;  
;  line 1: descriptor
;  line 2: crystal thickness (mm), crystal type
;  line 3: x points in scan, y points in scan, number of points in spectrum, crystal diameter (mm)
;  then there are nx by ny lines with the first 2 numbers being the x position and y position followed
;  by the spectrum 

	readf,1,desc                     ;read part description string
	readf,1,thickness,in_crystal     ;read part thickness (mm) and crystal type
	;  read scan points in x and y, number of points in spectrum and physical diameter of part (mm)
	readf,1,nx,ny,nspec,diam          

	in_crystal = strcompress(in_crystal,/remove_all)
	crystal = strlowcase(in_crystal)

 ;  set grating
 
  grating='1800'        ;grating used ('2400', '600', '1800', or 'echelle')
  if strpos(desc,'600') ge 0 then grating='600'
  
	if debug eq 'yes' then begin
		print,desc
		print,thickness
		print,'crystal:',crystal
		print,nx,ny,nspec
	endif


;  set up for plotting

  factor = fix(400./nx+0.5)   ;compute plotting scale factor
  factor = 10
  print,'factor:',factor
  
	bin = npixel/nspec		;CCD binning
	print,'bin:',bin

;   select dispersion for grating used

	case grating of

	'2400': disp = 0.00185d*bin 		   ;dispersion, (nm/pixel) 2400 l/mm grating
	'1800': disp = 0.00214d*bin       ;dispersion (nm/pixel) 1800 l/mm grating
	'600':  disp = 0.01144d*bin		   ;dispersion (nm/pixel) 600 l/mm grating
	'echelle': disp = 0.00119d*bin		 ;echelle dispersion, 632 nm in 9th order

	endcase

	wave = dindgen(nspec)*disp+wave0-double(nspec)*disp/2.d


;  nominal free spectral range and ld (nm)

	case crystal of
		'calcite': begin
						mu=calcite_biref(wave0,25.d)
						dlam = 10.d
						dmu = calcite_biref(wave0+dlam/2.d,25.d) - calcite_biref(wave0-dlam/2.d,25.d) 
						dmu_dlam = dmu/dlam         ;derivative of birefringence (nm-1)
						eff_mu=mu*(dmu_dlam*wave0/mu -1.d)			;effective birefringence (see Title, SoPh, 39, 505.)
						fsr=wave0^2/abs(eff_mu*thickness*1.d6)
						ld0=eff_mu*thickness*1.d6
					end

		'linb':	begin
						mu=linbo3_biref(wave0,25.d)
						mu2=linbo3_biref(wave0+1.d,25.d)
						dmu_dlam=mu2-mu     ;take derivative of birefringence (nm-1)
						eff_mu=mu*(dmu_dlam*wave0/mu -1.d)			;effective birefringence (see Title, SoPh, 39, 505.)
						fsr=wave0^2/abs(eff_mu*thickness*1.d6)
						ld0=eff_mu*thickness*1.d6
					end

		'quartz': begin
						mu=quartz_biref(wave0,25.d)
						mu2=quartz_biref(wave0+1.d,25.d)
						dmu_dlam=mu2-mu     ;take derivative of birefringence (nm-1)
						eff_mu=mu*(dmu_dlam*wave0/mu -1.d)			;effective birefringence (see Title, SoPh, 39, 505.)
						fsr=wave0^2/abs(eff_mu*thickness*1.d6)
						ld0=eff_mu*thickness*1.d6
					end

	endcase


	ifsr=fix(fsr/disp)			;convert free spectral range to pixels

	if debug eq 'yes' then print,fsr,ifsr

	ldelta=dblarr(nx,ny)
	intens=dblarr(nx,ny)
	x=dblarr(nx,ny)
	y=dblarr(nx,ny)


;  read in transmission spectra

  window,0,xs=690,ys=500,xpos=0,ypos=0
  window,6,xs=690,ys=500,xpos=0,ypos=0
  
  trans=dblarr(nspec)
	data=dblarr(nx,ny,nspec,/nozero)

	for iy=0,ny-1 do for ix=0,nx-1 do begin
		readf,1,xpos,ypos,trans
;		if debug eq 'yes' then print,ix,iy,xpos,ypos

		data[ix,ny-1-iy,*]=trans
		intens[ix,ny-1-iy]=max(trans)-min(trans)
		x[ix,iy]=xpos
		y[ix,iy]=ypos
	endfor
	close,1

	dx=abs(x[1,0]-x[0,0])
	dy=abs(y[0,1]-y[0,0])

	x=x-x[0,0]
	y=y-y[0,0]

;  fit transmission spectra

;  obtain guess for ldelta using grid search applied to spectrum near center of crystal

	num = 10000
	chisq = dblarr(num)

	int_thresh = 0.5			;intensity threshhold

	ld = ld0*( 0.85+0.3*(dindgen(num)/double(num-1)) )

	trans = data[nx/2,ny/2,*]			;get central spectrum
	wset,6
	plot,wave,trans

	rmin = run_min(trans,ifsr*1.)
	rmax = run_max(trans,ifsr*1.)
	trans = (trans-rmin)/(rmax-rmin)

	for j=0,num-1 do begin
		model = cos(!pi*ld[j]/wave)^2
		chisq[j] = total( (model-trans)^2 )
	endfor

	mn = min(chisq,imin)
	guess = [1.d0,ld[imin]]

	model = cos(!pi*ld[imin]/wave)^2
	wset,0
	plot,wave,trans,xsty=1,yr=[0,1],ysty=1,psym=3
	oplot,wave,model
	if debug eq 'yes' then read,'enter return',ans

	for iy=0,ny-1 do for ix=0,nx-1 do begin		;loop over all spectra

		if intens[ix,iy] gt int_thresh*max(intens) then begin
			trans = data[ix,iy,*]
			rmin = run_min(trans,ifsr*1.)
			rmax = run_max(trans,ifsr*1.)
			trans = (trans-rmin)/(rmax-rmin)

			if debug eq 'yes' then plot,wave,trans,xsty=1,yr=[0,1],ysty=1,psym=3


;  fit for length*delta with amoeba

			ftol = 1.d-8
			xfit = wave
			yfit = trans

			r = amoeba(ftol,function_name='fit_lyot',p0=guess,scale=10.,function_value=fval)

			model = r[0]*cos(!pi*r[1]/wave)^2
			if debug eq 'yes' then oplot,wave,model

			ldelta[ix,iy] = abs( r[1] )*1.d-6			;thickness * birefringence (mm)

			guess = [1.d0,r[1]]		;if not first time through, use last value as guess

		endif

	end


;  find center of sample

	bright = where(intens gt int_thresh*max(intens),complement=dark)		;identify bright data

	mask = dblarr(nx,ny)+1.d
	mask[dark] = 0.d
	xcent = total(x*mask)/total(mask)
	ycent = total(y*mask)/total(mask)


;  identify clear aperture region (fraction defined above as clear)

 	aperture = max(abs(x)) > max(abs(y))
  ;aperture = max(abs(x[bright])) > max(abs(y[bright]))
  r = sqrt( (x-xcent)^2+(y-ycent)^2 )
	clap_rad = clear*diam/2.d		;radius of clear aperture
	clap = where(r le clap_rad and intens gt int_thresh*max(intens),complement=bad)

	meanld = mean(ldelta[clap])
	minld = min(ldelta[clap])
	maxld = max(ldelta[clap])

	order = fix( 1.d6*ldelta/wave0 )  				;compute order across sample
	order0 = order[nx/2,ny/2]							;order at center of sample

	free = 1.d6*meanld*( 1.d/order0 - 1.d/(order0+1.d) ) 	;actual fsr at center (nm)
	print,'Free Spectral Range:',free,' nm'

	window,2,xs=factor*nx,ys=factor*ny,xpos=705,ypos=0,title='Intensity'
	tvscl,rebin(intens,factor*nx,factor*ny,/sample)

	window,3,xs=factor*nx,ys=factor*ny,xpos=720+factor*nx,ypos=0,title='Order'
	tvscl,rebin(order,factor*nx,factor*ny,/sample)


	tune = 1.d6*(ldelta/order)		     ;compute tuning wavelengths (nm)
	mean_tune = mean(tune[bright])     ;subtract mean to compute tuning wavelength variation (mean_tune is in nm)
	tune = tune-mean_tune
	tune[dark] = 0.d

  tune = tune/free          ;convert wavelength variation to units of the fsr
  
  
;  remove discontinuities defined as adjacent pixels differing by more than thresh waves

thresh=0.5d        

;  start at center and move in +/- x direction

	iy=ny/2
	for ix=nx/2,nx-1 do $
		if abs(tune[ix,iy]-tune[ix-1,iy]) gt thresh then $
			if tune[ix,iy] gt tune[ix-1,iy] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

	for ix=nx/2,0,-1 do $
		if abs(tune[ix,iy]-tune[ix+1,iy]) gt thresh then $
			if tune[ix,iy] gt tune[ix+1,iy] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

;  move in +/- y direction

	for ix=0,nx-1 do begin

		for iy=ny/2+1,ny-1 do $
		if abs(tune[ix,iy]-tune[ix,iy-1]) gt thresh then $
			if tune[ix,iy] gt tune[ix,iy-1] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

		for iy=ny/2-1,0,-1 do $
		if abs(tune[ix,iy]-tune[ix,iy+1]) gt thresh then $
			if tune[ix,iy] gt tune[ix,iy+1] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

	endfor

;  start at center and move in +/- y direction

	ix=nx/2
	for iy=ny/2,ny-1 do $
		if abs(tune[ix,iy]-tune[ix,iy-1]) gt thresh then $
			if tune[ix,iy] gt tune[ix,iy-1] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

	for iy=ny/2,0,-1 do $
		if abs(tune[ix,iy]-tune[ix,iy+1]) gt thresh then $
			if tune[ix,iy] gt tune[ix,iy+1] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

;  move in +/- x direction

	for iy=0,ny-1 do begin

		for ix=nx/2+1,nx-1 do if tune[ix,iy] ne 0. then $
		if abs(tune[ix,iy]-tune[ix-1,iy]) gt thresh then $
			if tune[ix,iy] gt tune[ix-1,iy] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

		for ix=nx/2-1,0,-1 do if tune[ix,iy] ne 0. then $
		if abs(tune[ix,iy]-tune[ix+1,iy]) gt thresh then $
			if tune[ix,iy] gt tune[ix+1,iy] then tune[ix,iy]=tune[ix,iy]-1 else tune[ix,iy]=tune[ix,iy]+1

	endfor

	tune=tune-mean(tune[bright])
	tune[dark]=0.d


;  optionally remove low order fit from wavelength shift

	xfit=x-x[nx/2,ny/2]
	yfit=y-y[nx/2,ny/2]
	mx=max(abs(x)) > max(abs(y))
	xfit=xfit/mx
	yfit=yfit/mx

	if rem_zernike eq 'yes' then begin
		coefs=zernike_fit(xfit[bright],yfit[bright],tune[bright],nz,tune_fit)
		tune[bright]=tune[bright]-tune_fit
		tune=tune-median(tune[bright])
		tune[dark]=0.d
	endif

	if rem_poly eq 'yes' then begin
		coefs=surf_fit_xy(tune[bright],xfit[bright],yfit[bright],np)
		fit=eval_surf(coefs,xfit[bright],yfit[bright])
		tune[bright]=tune[bright]-fit
		tune=tune-median(tune[bright])
		tune[dark]=0.d
	endif

	if rem_trend eq 'yes' then begin
		coefs=surf_fit_xy(tune[bright],xfit[bright],yfit[bright],np)
		fit=eval_surf(coefs,xfit,yfit)
		tune[dark]=fit[dark]		;fill outside with polynomial fit
		;sm=smooth(tune,nx/5,/edge_truncate)
		sm=smooth(tune,5,/edge_truncate)
		tune=tune-sm
		tune=tune-median(tune[bright])
		tune[dark]=0.d
	endif


;  compute rms and ptv variations across crystal clear aperture

	rms=stdev(tune[clap])							;rms variation over clear aperture (units of fsr)
	ptv=max(tune[clap])-min(tune[clap])		;peak-to-valley over clear aperture (units of fsr)


;  compute histogram

	window,4,xs=400,ys=400,xpos=1540
	h=histogram(tune[clap],min=min(tune[clap]),max=max(tune[clap]),nbins=20,locations=xvals)
	plot,xvals,h,psym=10,xtit='Wavelength Shift (waves)',ytit='Frequency',chars=1

	fw=fwhm(h,xvals)		;fwhm of distribution (units of fsr)


;  display variation of bandpass shift in waves

	rtune=rebin(tune,factor*nx,factor*ny,/sample)

	display_crystal_image,rtune,-0.7*ptv,0.7*ptv,tit='Wavelength Shift', $
	 windnum=1,table=3,output='s',wherebar='right',units='(waves @ 635nm)', thk=2, $
	 xt='X Position (mm)',yt='Y Position (mm)',xscale=dx/factor,yscale=dy/factor,offset=0.

;  draw clap

	angle=2.*!pi*dindgen(100)/99.
	xline=clap_rad*sin(angle)+xcent+0.5
	yline=clap_rad*cos(angle)+ycent+0.5
	oplot,xline,yline,color=255,thick=2


;  annotate plot

	wset,1
  thick,2
	date_time=strmid(files[ifile],18,15,/reverse_offset)
	outfile=date_time+' '+desc+'_waves.bmp' 

	xyouts,0.1,0.92,/norm,desc,color=0,chars=1.5
	xyouts,0.1,0.88,/norm,date_time,color=0,chars=1.5
	str=in_crystal+string(format='(4x,f6.2," mm Thick")',thickness)
	xyouts,0.1,0.84,/norm,str,color=0,chars=1.5
	str=string(format='("Measured FSR:",f6.3," nm")',free)
	xyouts,0.1,0.80,/norm,str,color=0,chars=1.5

	str=string(format='("CLAP:",f6.1," mm")',clap_rad*2.)
	xyouts,0.1,0.18,/norm,str,color=0,chars=1.5
  str=string(format='("RMS:",f6.3," nm,",f6.3," FSR")',rms*free,rms)
	xyouts,0.1,0.14,/norm,str,color=0,chars=1.5
	str=string(format='("PTV:",f6.3," nm,",f6.3," FSR")',ptv*free,ptv)
	xyouts,0.1,0.10,/norm,str,color=0,chars=1.5
	str=string(format='("FWHM:",f6.3," nm,",f6.3," FSR")',fw*free,fw)
	xyouts,0.1,0.06,/norm,str,color=0,chars=1.5

	if rem_poly eq 'yes' then begin
		str=string(format='("Polynomial Order:",i2)',np)
		xyouts,0.1,0.02,/norm,str,color=0,chars=1.5
		outfile=date_time+' '+desc+' -P.bmp'
	endif

	if rem_zernike eq 'yes' then begin
		str=string(format='("Zernike Polynomials:",i2)',nz)
		xyouts,0.1,0.02,/norm,str,color=0,chars=1.5
		outfile=date_time+' '+desc+' -Z.bmp'
	endif

	if rem_trend eq 'yes' then begin
		str='Trend Removed'
		xyouts,0.1,0.02,/norm,str,color=0,chars=1.5
	endif

;  write bmp

  outfile=date_time+' '+desc+'_waves.bmp'
  if rem_trend eq 'yes' then outfile=date_time+' '+desc+'_waves-T.bmp'

	write_bmp,out_dir+outfile,tvrd(true=1),r,g,b,/rgb


;  compute birefringence variation from equation dmu/mu = 2 dlambda/lambda
;  must convert wavelength variation (tune) from waves to wavelength units (multiply by free=fsr)

  dmu=2.d*abs(mu)*free*tune/mean_tune
  dmu=dmu-mean(dmu[bright])
  dmu[dark]=0.d

;  save dmu map to idl savefile

  savfile=date_time+' '+desc+'_biref.sav'  
  if rem_trend eq 'yes' then savfile=date_time+' '+desc+'_biref-T.sav'
  save,dmu,tune,file=out_dir+savfile
  
;  plot birefringence variation across crystal

  reb_dmu=rebin(dmu,factor*nx,factor*ny,/sample)
  
  rms=stdev(dmu[clap])             ;rms variation over clear aperture
  
  ptv=max(dmu[clap])-min(dmu[clap])   ;peak-to-valley over clear aperture
  
  display_crystal_image,reb_dmu*1.e5,-0.7e5*ptv,0.7e5*ptv,tit='Birefringence Map', $
    windnum=7,table=3,output='s',wherebar='right',units='10!e-5!n', thk=2, $
    xt='X Position (mm)',yt='Y Position (mm)',xscale=dx/factor,yscale=dy/factor,offset=0.


;display_crystal_image_2,reb_dmu*1.e5,-1.,1.,tit='Birefringence Map', $
;  windnum=7,table=3,output='s',wherebar='right',units='x10!e-5!n', thk=4, $
;  xt='X Position (mm)',yt='Y Position (mm)',xscale=dx/factor,yscale=dy/factor,offset=0.

;  draw clap
  
  angle=2.*!pi*dindgen(100)/99.d
  xline=clap_rad*sin(angle)+xcent+0.5d
  yline=clap_rad*cos(angle)+ycent+0.5d
  oplot,xline,yline,color=255,thick=2
  
;  annotate plot
  
  wset,7
  thick,2
  date_time=strmid(files[ifile],18,15,/reverse_offset)

  do_anno='yes'          ;annotate? ('yes' or 'no')
  if do_anno eq 'yes' then begin
  
  xyouts,0.1,0.92,/norm,desc,color=0,chars=1.5
  xyouts,0.1,0.88,/norm,date_time,color=0,chars=1.5
  str=in_crystal+string(format='(4x,f6.2," mm Thick")',thickness)
  xyouts,0.1,0.84,/norm,str,color=0,chars=1.5
  str=string(format='("Measured FSR:",f6.3," nm")',free)
  xyouts,0.1,0.80,/norm,str,color=0,chars=1.5
  
  str=string(format='("CLAP:",f6.1," mm")',clap_rad*2.)
  xyouts,0.1,0.18,/norm,str,color=0,chars=1.5
  str=string(format='("RMS:",f7.2," *10!e-6!n")',rms*1.e6)
  xyouts,0.1,0.14,/norm,str,color=0,chars=1.5
  str=string(format='("PTV:",f7.2," *10!e-6!n")',ptv*1.e6)
  xyouts,0.1,0.10,/norm,str,color=0,chars=1.5
  str=string(format='("FWHM:",f6.3," nm,",f6.3," FSR")',fw*free,fw)
  ;xyouts,0.1,0.08,/norm,str,color=0,chars=1.5

  if rem_trend eq 'yes' then begin
    str='Trend Removed'
    xyouts,0.1,0.02,/norm,str,color=0,chars=1.5
  endif

  thick,0
  endif
  
;  write bmp
  
  outfile=date_time+' '+desc+'_biref.bmp'
  if rem_trend eq 'yes' then outfile=date_time+' '+desc+'_biref-T.bmp'

  write_bmp,out_dir+outfile,tvrd(true=1),r,g,b,/rgb


;  compute transmission profile averaged over clear aperture

  n=200         ;number of wavelengths
  dlam=1.4d      ;wavelength range (nm)
  w=dlam*dindgen(n)/double(n-1)-dlam/2.+1075.7d
  dlam=0.3d       ;wavelength range (nm)
  w=dlam*dindgen(n)/double(n-1)-dlam/2.+530.36d
  dlam=0.8d       ;wavelength range (nm)
  w=dlam*dindgen(n)/double(n-1)-dlam/2.+657.08d

  avg_trans=dblarr(n)
  avg_trans_uni=dblarr(n)
  
  th=[44.,22.,11.,5.5,2.75]*1.d6   ;stage thicknesses (nm)
  ns=n_elements(th)
  dmu=dmu-mean(dmu[clap])
  
  for k=0,n-1 do begin 
    total_mu=dmu+mu                 ;compute measured birefringence from fractional birefringence
    trans=dblarr(nx,ny)+1.d          ;transmission with non-uniform birefringence
    trans_uni=dblarr(nx,ny)+1.d      ;transmission using uniform birefringence
    for i=0,ns-1 do begin
      trans=trans*cos(!pi*total_mu*th[i]/w[k])^2   ;transmission across element
      trans_uni=trans_uni*cos(!pi*mu*th[i]/w[k])^2   ;theoretical transmission across element
    endfor
    avg_trans[k]=mean(trans[clap])
    avg_trans_uni[k]=mean(trans_uni[clap])
  endfor
  print,fwhm(avg_trans_uni,w),fwhm(avg_trans,w)
  
  window,5,xs=600,ys=400,xpos=1950,ypos=0
  plot,w,avg_trans_uni,xtit='Wavelength (nm)',ytit='Transmission',chars=1.5,xsty=1
  oplot,w,avg_trans

	wait,0.005
endfor

print,'done'
end

function fit_lyot,p

common fit, xfit,yfit

model=p[0]*cos(!pi*p[1]/xfit)^2
chisq=total( (model-yfit)^2 )
return,chisq
end

function ldfunc, params

common fit, xfit,yfit

x=( xfit-min(xfit) )/( max(xfit)-min(xfit) )
f=params[0] + params[1]*x +	params[2]/(x-params[3]) + params[4]*x^3

return,f
end

function ldmin, params

common fit, xfit,yfit

return,total( (call_function('ldfunc',params)-yfit)^2 )

end
