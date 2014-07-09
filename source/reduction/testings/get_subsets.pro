pro get_subsets

; make subset of images (at different temperatures) for testing 
; 
Compile_opt idl2

!p.background=cgcolor('white')
plotsym, 0, /fill, 0.9
window, 0, xsize=600, ysize=550, xpos=1000, ypos=100

indir='~/BRITE/data/UB/p1/ORION/'

outdir='~/BRITE/data/UB/testing/ORION/p1_subset1/'

filesin=file_search(indir+'*.sav', count=nf)

for i=0, nf-1 do begin
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                      ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype
                      
                      
                      jd1=jd-jd[0]
                      
                      plot, jd1, roi_dim[0,*], color=cgcolor('black'), psym=2, /ynozero
                      
                      
                      stop
                      
  jd1=jd-jd[0]
                      
  xx=where(finite(ccd_temp) eq 1, nxx)
  
  ccd_temp=ccd_temp[xx]
  jd1=jd1[xx]
  medimg0=medimg0[xx]
  data1=data1[*,*,xx]
    
  cghistoplot, ccd_temp, binsize=1, histdata=result, locations=loc, mininput=10, reverse_indices=ri
  
  n=10
  ntemp=n_elements(loc)
  
  if i eq 0 then begin
    
  
  for j=0, ntemp-1 do begin
    
    iloc=round(randomu(seed, n)*result[j])
    
    allri=ri[ri[j]:ri[j+1]-1]
    
    if j eq 0 then sub1=allri[iloc] else sub1=[sub1, allri[iloc]]
    
  endfor
  
  endif else sub1=sub1
  
  ccd_temp1=ccd_temp[sub1]
  jd2=jd1[sub1]
  medimg=medimg0[sub1]
  data2=data1[*,*,sub1]
  jd=jd[xx[sub1]]
  
  data1=data2
  ccd_temp=ccd_temp1
  medimg0=medimg
   
  ; save new file
  fileout=outdir+roi_name[0]+'_ss1.sav'
  
  save, filename=fileout, jd, roi_name, data1, ccd_temp, medimg, xc, yc, vmag
  
endfor

print, 'end of program'

end

