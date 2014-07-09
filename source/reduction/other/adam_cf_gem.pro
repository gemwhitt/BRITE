pro adam_cf_gem

; compare Adam LCs with Gemma LCs
; 
Compile_opt idl2

adam_dir='/Users/gemmawhittaker/BRITE/data/UB/adam_lcs/'
gem_dir='/Users/gemmawhittaker/BRITE/data/UB/p4_2/'

afiles=file_search(adam_dir+'HD*.txt', count=nafl)

gfiles=file_search(gem_dir+'HD*.sav', count=ngfl)

for i=0, nafl-1 do begin
  
  readcol, afiles[i], ta, fa, ea, format='d,d,d'
  
  restore, gfiles[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
  tg=jd
  fg=(-2.5)*alog10(flux)
  
  npta=n_elements(ta)
  nptg=n_elements(jd)
  
  print, file_basename(afiles[i], '.txt'), ', ', file_basename(gfiles[i], '.sav')
  print, npta, nptg
  print, robust_mean(fa,3), robust_mean(fg,3)
  print, ''
  
  stop
  
  
  
endfor

stop
print, 'End of program'

end