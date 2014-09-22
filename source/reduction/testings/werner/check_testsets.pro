pro check_testsets

Compile_opt idl2

indir='/Users/gemmawhittaker/BRITE/TESTSETS/werner4lc/p2/'
;indir='/Users/gemmawhittaker/BRITE/UB/CENTAURUS/data/p2/'

filein=file_search(indir+'HD127973*TO*.sav', count=nfile)

for f=0, nfile-1 do begin
  
  fname=file_basename(filein[f], '.sav')
  print, fname
  
  restore, filein[f]  ;jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag
  
  nimg=n_elements(jd)
  print, nimg, ' number of images'
  
  ; date range
  caldat, jd[0], mon, day, yr
  print, 'Start date is '+strtrim(day,2)+' '+strtrim(mon,2)+' '+strtrim(yr,2)
  caldat, jd[nimg-1], mon, day, yr
  print, 'End date is '+strtrim(day,2)+' '+strtrim(mon,2)+' '+strtrim(yr,2)
  
  ;plotsym, 0, /fill, 0.8
  ;plot, jd, medimg0, color=cgcolor('black'), psym=8
  
  ;xx=where(flag eq 0, nxx)
  ;if nxx gt 0 then oplot, jd[xx], medimg0[xx], color=cgcolor('green'), psym=8
  
  for im=1020, 1050 do begin  ;nimg-1 do begin
    
    print, im
    plot_image, bytscl(data1[*,*,im], 20, 500), title=strtrim(ccd_temp[0,im],2)+', '+strtrim(jd[im],2),$
     color=cgcolor('black'), charsize=1.
    stop
    
  endfor
  stop
  
endfor

print, 'Re-do medcol removals? OR go straight to photometry'

end