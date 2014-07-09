pro check_nobin_results

; program written to check the results from p3_p4_custom_4 - using no binning
; 
; check mags, number of counts, flux, residual, tracking...
;
Compile_opt idl2

!p.background=cgcolor('white')

indir='/Users/gemmawhittaker/BRITE/data/UB/p4_nobin/'

filesin=file_search(indir+'HD35715_nobin.sav', count=nfiles)

if nfiles eq 0 then stop

window, 0, xsize=600, ysize=350, xpos=100, ypos=400
window, 1, xsize=600, ysize=350, xpos=800, ypos=400
window, 2, xsize=600, ysize=350, xpos=300, ypos=0

plotsym, 0, 0.7, /fill

dcount=fltarr(nfiles)
fscat=fltarr(nfiles)
rscat=fltarr(nfiles)
mag=fltarr(nfiles)


for i=0, nfiles-1 do begin
  
  restore, filesin[i] ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, dresid, counts, vmag
  
  mag[i]=vmag
  
  xx=where(jd gt 0.02 and jd lt 0.09, nxx)
  
  flux1=flux[xx]
  time1=jd[xx]
  count=counts[xx]
  resid=dresid[xx]
  
  print, file_basename(filesin[i], '.sav')
  print, vmag
  
  ;flux1=flux/robust_mean(flux,3)
  
  wset, 0
  plot, time1, flux1, color=cgcolor('black'), psym=8, xtitle='Time', ytitle='Flux', /ynozero, $
    charsize=0.9;, xrange=[0.56,0.57]
    
    wset, 1
  plot, time1, resid, color=cgcolor('black'), psym=8, xtitle='Time', ytitle='Resid', /ynozero, $
    charsize=0.9;, xrange=[0.56,0.57]
    
    wset, 2
  plot, time1, count, color=cgcolor('black'), psym=8, xtitle='Time', ytitle='Counts', /ynozero, $
    charsize=0.9;, xrange=[0.56,0.57]
    
    dcount[i]=max(count)-min(count)
    fscat[i]=stddev(flux1/robust_mean(flux1,3))
    rscat[i]=stddev(resid/robust_mean(resid,3))
  
  stop
endfor
wset, 0
plot, mag, dcount, color=cgcolor('black'), psym=8, xtitle='Magnitude', ytitle='Delta Count'

wset, 1
plot, mag, rscat, color=cgcolor('black'), psym=8, xtitle='Magnitude', ytitle='Resid STDDEV'

;stop
wset, 2
plot, mag, fscat, color=cgcolor('black'), psym=8, xtitle='Magnitude', ytitle='Flux STDDEV'



stop
print, 'End of program'
end

