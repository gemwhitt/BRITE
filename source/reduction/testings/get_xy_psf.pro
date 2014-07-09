pro get_xy_psf

Compile_opt idl2

; compare slavek PSF-fitting xy-locations with brite_findpeaks locations....

indir1='~/BRITE/data/UB/p4_sr1/'

indir2='~/BRITE/data/UB/p4_gem1/sav/'

filesin1=file_search(indir1+'*', count=nf1)
filesin2=file_search(indir2+'Orion-CF1-2*', count=nf2)

for i=0, nf1-1 do begin
  
  restore, filesin1[i]
  print, file_basename(filesin1[i],'.sav')
  
  xsl=reform(xy_psf[0,*])
  ysl=reform(xy_psf[1,*])
  
  stop
  
  
  restore, filesin2[i]
  print, file_basename(filesin2[i],'.sav')
  
  xg=reform(xy_psf[0,*])
  yg=reform(xy_psf[1,*])
  
  flux=flux[0,*]
  
 ; plot, xsl, 

  
  stop
  
  
  
endfor

stop




end