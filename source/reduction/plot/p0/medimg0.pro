pro medimg0

; PURPOSE: Investigate medimg0 - calculated in save_medcol
; Output - make plots - time and medimg0, ccd_temp and medimg0
; 
; 
Compile_opt idl2
plotsym, 0, /fill, 0.3

makeplot='y'

sat='UB'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'

outdir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/plots/p0/'

statfile='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/stats/p0/ccd_loc_vs_medimg0.txt'

tardirs=file_search(indir+'HD*/', count=ntar) ; array of target directories containing .sav files

for tar=0, ntar-1 do begin
  
  filesin=file_search(tardirs[tar]+'/*.sav', count=nf)
  
  fname=file_basename(tardirs[tar])
  
  ; make arrays for storing variables for the whole observation
  time1=[]
  medimg=[]
  temps=[]
  xloc=[]
  yloc=[]
  
  for f=0, nf-1 do begin
    
    obj=obj_new('IDL_Savefile', filesin[f])
    obj->restore, 'jd'
    obj->restore, 'medimg0'
    obj->restore, 'ccd_temp'
    obj->restore, 'roi_loc'
    obj->restore, 'exp_ttl'
    obj->restore, 'exp_time'
    
    ;only save info if 1s non-stacked exposures
    if exp_ttl[0] ne 1. OR exp_time[0] ne 1. then continue
    
    time1=[time1,jd]
    
    nimg=n_elements(jd)
    
    ; calculate average CCD_temp
    avg_temp=fltarr(nimg)
    for im=0, nimg-1 do avg_temp[im]=average(ccd_temp[*,im])
    temps=[temps,avg_temp]
    
    ; calculate x and y CCD pixel
    xloc=[xloc,reform(roi_loc[0,*])]
    yloc=[yloc,reform(roi_loc[2,*])]
    
    ; calculate the median medimg0 and average medimg0
    medimg=[medimg,medimg0]
    
  endfor
  
  ; sort by time increasing
  sort1=sort(time1)
  time1=time1[sort1]
  medimg=medimg[sort1]
  temps=temps[sort1]
  xloc=xloc[sort1]
  yloc=yloc[sort1]
  
  jd1=time1-time1[0]
  
  caldat, time1[0], mon, day, yr
  startdate=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
  
  if makeplot eq 'y' then begin
    plotout=outdir+fname+'_medimg0.ps'
    ps_on, plotout, xsize=16, ysize=26
    !p.multi=[0,1,2,0,0]
    plot, jd1, medimg, color=cgcolor('black'), psym=8, xtitle='Time in days since '+startdate, $
      ytitle='Median of raw image', title=fname, charsize=0.8
    plot, jd1, medimg, color=cgcolor('black'), psym=8, xtitle='Time in days since '+startdate, $
      ytitle='Median of raw image', yrange=[0,200], charsize=0.8     
    plot, avg_temp, medimg, color=cgcolor('black'), psym=8, xtitle='Average CCD temperature', $
      ytitle='Median of raw image', title=fname, charsize=0.8
    plot, avg_temp, medimg, color=cgcolor('black'), psym=8, xtitle='Average CCD temperature', $
      ytitle='Median of raw image', yrange=[0,200], charsize=0.8      
    ps_off
    !p.multi=0
    !p.background=cgcolor('white')
  endif
  
  xloc=median(xloc)
  yloc=median(yloc)
  med_img0=median(medimg)
  
  ; write out read_dist, med_img0 and avg_img0 to file fo all targets
  openw, lun, statfile, /get_lun, /append
  printf, lun, fname, xloc, yloc, med_img0, format='(a,x,i,x,i,x,i)'
  free_lun, lun
  
  
endfor

print, 'end of program'
end