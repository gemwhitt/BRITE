pro p0_plots

; PURPOSE: plot medimg0 and ccd_temp vs time 
; 
Compile_opt idl2
plotsym, 0, /fill, 1.1

sat='UB'
field='CENTAURUS'

indir='~/BRITE/data/'+sat+'/p1/'+field+'/'

outdir='~/BRITE/data/'+sat+'/reduction/'+field+'/p0_plots/'

filesin=file_search(indir+'*.sav', count=nf)

; calculate median of medimg0's and plot versus x and y on CCD
med_medimg=lonarr(nf)
xccd=fltarr(nf)
yccd=fltarr(nf)

for fl=0, nf-1 do begin
  
  obj=obj_new('IDL_Savefile', filesin[fl])
  obj->restore, 'jd'
  ;obj->restore, 'medimg0'
  obj->restore, 'medimg'
  obj->restore, 'ccd_temp'
  ;obj->restore, 'roi_loc'
  obj->restore, 'roi_dim'
  
  roi_loc=roi_dim
  medimg0=medimg
  
  fname=file_basename(filesin[fl], '_p1.sav')
  
  jd1=jd-jd[0]
  nfrm=n_elements(jd1)
  
  caldat, jd[0], mon, day ,yr
  date1=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
  
  ; calculate the duty cycle
  totdur=jd1[nfrm-1]

  ; get cadence of single and stacked images
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  cadence=median(jdiff)                                                     
  cadence_sec=cadence*24.*60.*60                                                ; SAVE - in secs
     
  ; get percentage of time observed......
  ; calculate expected number of points based on exposure time and 15 min orbit every 100 mins
  tot_obs_time=(nfrm*cadence_sec/60.)  ; mins
  possible_obs_time=totdur*24.*60./100.*15.                 ; mins
  duty_cyc=fix(tot_obs_time/possible_obs_time*100.)              ; SAVE - % number of points observed in totdur
  
  ; determine number of "saturated" frames i.e with medimg0 gt 5000
  bad=where(medimg0 ge 5000, nbad, complement=good)
  ; convert this to a %
  badpc=round(float(nbad)/float(nfrm)*100.)
  
  med_medimg[fl]=median(medimg0[good])
  xccd[fl]=roi_loc[0]
  yccd[fl]=roi_loc[2]
  
  ;plotout=outdir+fname+'_medimg0.ps' 
  ;ps_on, plotout, xsize=18, ysize=26
  ;!p.multi=[0,1,2,0,0]
  ;plot, jd1, medimg0, color=cgcolor('black'), psym=8, xtitle='Days since '+date1, ytitle='p0 Image median', $
  ; title=fname, charsize=0.8
  ;oplot, [0, jd1[nfrm-1]+10.], [16400,16400], thick=3, color=cgcolor('red')
  ;oplot, [0, jd1[nfrm-1]+10.], [5000,5000], thick=3, color=cgcolor('purple')
  ;xyouts, 1, 2000, 'Duty cycle = '+strtrim(duty_cyc,2)+'%', color=cgcolor('blue'), charsize=0.7
  ;xyouts, 1, 15000, 'Bad frms = '+strtrim(badpc,2)+'%', color=cgcolor('red'), charsize=0.7
  ;plot, jd1, medimg0, color=cgcolor('black'), psym=8, xtitle='Days since '+date1, ytitle='p0 Image median', $
  ;  title=fname, charsize=0.8, yrange=[50, 150]
  ;ps_off
  
  ;if fl eq 0 then begin
  ;  plotout=outdir+'time_temp.ps'
  ;  ps_on, plotout, xsize=18, ysize=12
  ;  !p.multi=0
  ;  plot, jd1, ccd_temp[0,*], color=cgcolor('black'), psym=8, xtitle='Days since '+date1, ytitle='CCD Temperature', $
  ;    title='ORION BA', charsize=0.8
  ;  ps_off
  ;endif
  
endfor
wset, 0
plot, xccd, med_medimg, color=cgcolor('black'), psym=8, yrange=[90,190], xrange=[0,4000], xstyle=1, $
  xtitle='X-position of raster on CCD', ytitle='Average median image', charsize=0.9
oplot, [0, 4000], [median(med_medimg),median(med_medimg)], thick=3, color=cgcolor('blue')
oplot, [2000, 2000], [0,200], thick=3, color=cgcolor('purple'), linestyle=2

wset, 1
plot, yccd, med_medimg, color=cgcolor('black'), psym=8, yrange=[90,190], xrange=[0,2700], xstyle=1, $
  xtitle='Y-position of raster on CCD', ytitle='Average median image', charsize=0.9
oplot, [0, 4000], [median(med_medimg),median(med_medimg)], thick=3, color=cgcolor('blue')
oplot, [1350, 1350], [0,200], thick=3, color=cgcolor('purple'), linestyle=2
stop
print, 'end of program'
end