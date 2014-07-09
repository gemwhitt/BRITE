pro compare_results_bins

; program to compare the results from PSF-fitting, using the p3 images (HP and orbit-trend removed)
; fitting was tried at 4x4,8x8 and 16x16 binning - computation time was <1, 3, 13 mins, respectively.
; 
Compile_opt idl2

!p.background=cgcolor('white')

indir='~/BRITE/data/UB/testing/nostack/p4/'

filesin=file_search(indir+['HD34085_rb8_']+'*.sav', count=nfiles)

fname=file_basename(filesin, '.sav')

pos1=strpos(fname, '_')

hdname=strarr(nfiles)
for i=0, nfiles-1 do hdname[i]=strmid(fname[i], 0, pos1[i])
targets=hdname[uniq(hdname, sort(hdname))]

for i=0, n_elements(targets)-1 do begin
  
  files1=filesin[where(hdname eq targets[i], nfiles)]
  
  obj1=obj_new('IDL_Savefile', files1[0])
  
  obj1->Restore, 'jd' 
  
 ; npts=n_elements(jd)
  npts=516
  
  fl=fltarr(nfiles, npts)
  fle=fltarr(nfiles, npts)
  bkg=fltarr(nfiles, npts)
  bke=fltarr(nfiles, npts)
  
  for j=0, nfiles-1 do begin
    
    restore, files1[j]
    
    fl[j,*]=flux[0:515]
    fle[j,*]=ferr[0:515]
    bkg[j,*]=bkgd[0:515]
    bke[j,*]=berr[0:515]
    
  endfor

; compute stats
mean_lc=robust_mean(fl,2)
scat_lc=robust_sigma(fl)
mean_bk=robust_mean(bkg,2)
scat_bk=robust_sigma(bke)  

jd1=jd-jd[0]
; find gaps in orbit
jd2=jd1[1:npts-1]
jdiff=jd2-jd1
gap=where(jdiff gt 0.015, ngap)
gap=[gap,npts-1]

; do plots
cols1=cgcolor(['blue', 'green', 'red'])

window, 0, xsize=700, ysize=500

for kk=0, ngap-1 do begin
  
  if kk eq 0 then begin
    time=jd1[0:gap[0]]
    flux=fl[*,0:gap[0]]
  endif else begin
    time=jd1[gap[kk]+1:gap[kk+1]]
    flux=fl[*,gap[kk]+1:gap[kk+1]]
  endelse
  
  plot, time, flux[0,*], color=cgcolor('black'), xtitle='Time (days)', ytitle='Flux', charsize=0.7, $
    title=targets[i], yrange=[0.98,1.02]
  for j=0, nfiles-1 do oplot, time, flux[j,*], color=cols1[j]
    wait, 2
endfor

stop
endfor

stop
end