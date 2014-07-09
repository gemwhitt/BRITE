pro make_ascii

Compile_opt idl2
  
indir='/Users/gemmawhittaker/BRITE/data/UB/p4_3/'
outdir='/Users/gemmawhittaker/BRITE/data/UB/p4_3/'
  
filesin=file_search(indir+'*.sav', count=nf)
  
for i=0, nf-1 do begin
  
  restore, filesin[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
  n=n_elements(jd)
  
  fname=file_basename(filesin[i], '_3.sav')
  
;  stop
  
  fileout=outdir+fname+'_v2.txt'
  
  xloc=xy_psf[0,*]+roi_dim[0]
  yloc=xy_psf[1,*]+roi_dim[2]
  
  mags=-2.5*alog10(flux)
  
  openw, lun, fileout, /get_lun
  for j=0, n-1 do printf, lun, jd[j], flux[j], mags[j], xloc[j], yloc[j], ccd_temp[j], format='(d,x,f,x,f,x,f,x,f,x,f)'
  free_lun, lun
  
endfor

print, 'end of program'

end
  