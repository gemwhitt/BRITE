pro orion_ba

Compile_opt idl2

indir='/Users/gemmawhittaker/BRITE/data/BA/ORION/p1/sav/'

plotout='/Users/gemmawhittaker/BRITE/data/BA/ORION/plots/ba_orion_ccdtemp.ps'
pdf='/Users/gemmawhittaker/BRITE/data/BA/ORION/plots/ba_orion_ccdtemp.pdf'

txtout='/Users/gemmawhittaker/BRITE/results/ba_ub/ORION/ba_temps.txt'

filesin=file_search(indir+'HD31237*.sav', count=nf)

jd0=[]
temps=[]
ii=[]

for ff=0, nf-1 do begin
  
  restore, filesin[ff]  ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
                        ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, data1
                        
  nimg=n_elements(jd)
  
  ii=[ii,(intarr(nimg)+ff)] ; ii will be 0 or 1
                        
  jd0=[jd0,jd]    
  
  ccd_temp2=fltarr(nimg)
  for im=0, nimg-1 do ccd_temp2[im]=average(ccd_temp[*,im])  
  
  temps=[temps,ccd_temp2]                 
  
endfor

; sort arrays - increasing JD
sort1=sort(jd0)

jd0=jd0[sort1]
temps=temps[sort1]
ii=ii[sort1]

; print out these variables in .txt file
nimg=n_elements(jd0)

openw, lun, txtout, /get_lun
for ii=0, nimg-1 do printf, lun, jd0[ii], temps[ii], format='(d14.6,x,f7.3)'
free_lun, lun

jd1=jd0-jd0[0]

; get start date of observations
caldat, jd0[0], mon, day, yr
stdate=strtrim(day,2)+'/'+strtrim(mon,2)+'/'+strtrim(yr,2)

;plotsym, 0, /fill, 0.2
;ps_on, plotout, xsize=15, ysize=11
;plot, jd1, temps, color=cgcolor('black'), title='BA ORION', xtitle='Days from '+stdate, ytitle='Average CCD temp', $
;  charsize=0.7, /nodata
;oplot, jd1[where(ii eq 1)], temps[where(ii eq 1)], psym=8, color=cgcolor('sea green')
;oplot, jd1[where(ii eq 0)], temps[where(ii eq 0)], psym=8, color=cgcolor('purple')
;ps_off

;spawn, 'convert '+plotout+' '+pdf


stop
end

