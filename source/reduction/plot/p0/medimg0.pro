pro medimg0

; PURPOSE: Investigate medimg0 - calculated in save_medcol
; Output - make plots - time and medimg0, ccd_temp and medimg0
; 
; 
Compile_opt idl2
plotsym, 0, /fill, 0.3

sat='BA'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'

outdir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/plots/p0_images/'

filesin=file_search(indir+'*p0.sav', count=nf)

statfile='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/stats/p0_images/ccd_loc_vs_medimg0.txt'

for f=0, nf-1 do begin
  
  fname=file_basename(filesin[f], '_p0.sav')
  
  obj=obj_new('IDL_Savefile', filesin[f])
  obj->restore, 'jd'
  obj->restore, 'medimg0'
  obj->restore, 'ccd_temp'
  obj->restore, 'roi_loc'
    
  jd1=jd-jd[0]
  
  nimg=n_elements(jd1)
  
  caldat, jd, mon, day ,yr
  
  startdate=strtrim(day[0],2)+'_'+strtrim(mon[0],2)+'_'+strtrim(yr[0],2)
  
  ; calculate average CCD_temp
  avg_temp=fltarr(nimg)
  for im=0, nimg-1 do avg_temp[im]=average(ccd_temp[*,im])
  
  ; calculate distance from readout using roi_loc
  read_dist=roi_loc[0,0]+roi_loc[2,0]
  
  ; calculate the median medimg0 and average medimg0
  med_img0=median(medimg0)
  avg_img0=average(medimg0)
  
  ; write out read_dist, med_img0 and avg_img0 to file fo all targets
  openw, lun, statfile, /get_lun, /append
  printf, lun, fname, read_dist, roi_loc[0,0], roi_loc[2,0], med_img0, avg_img0, format='(a,x,i,x,i,x,i,x,i,x,f)'
  free_lun, lun
  stop
  continue
  
  fout1=outdir+fname+'_medimg0.ps'
  ps_on, fout1, xsize=16, ysize=26
  !p.multi=[0,1,2,0,0]
  plot, jd1, medimg0, color=cgcolor('black'), psym=8, xtitle='Time in days since '+startdate, $
    ytitle='Median of raw image', title=fname, charsize=0.8
  plot, jd1, medimg0, color=cgcolor('black'), psym=8, xtitle='Time in days since '+startdate, $
    ytitle='Median of raw image', yrange=[0,200], charsize=0.8
  
  plot, avg_temp, medimg0, color=cgcolor('black'), psym=8, xtitle='Average CCD temperature', $
    ytitle='Median of raw image', title=fname, charsize=0.8
  plot, avg_temp, medimg0, color=cgcolor('black'), psym=8, xtitle='Average CCD temperature', $
    ytitle='Median of raw image', yrange=[0,200], charsize=0.8
    
  ps_off
  
endfor


print, 'end of program'
end