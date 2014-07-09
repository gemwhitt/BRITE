pro p1_images

; check p1 images try to establish trends and thresholds for image categories...
;
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 0.9
;window, 0, xsize=600, ysize=550, xpos=1000, ypos=200
 
; use test data
;indir='~/BRITE/data/UB/testing/ORION/p1_subset1/'
indir='~/BRITE/data/UB/p1/ORION/'

outdir='~/BRITE/data/UB/plots/ORION/systematic_trends/'
  
filesin=file_search(indir+'*.sav', count=nf)

for i=0, nf-1 do begin
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                      ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype
                      
  print, roi_name
  
  jd1=jd-jd[0]
  
  plot, jd1, ccd_temp, color=cgcolor('black'), psym=2
  stop

  nfrm=n_elements(jd)
  
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  totdn=lonarr(nfrm)
  
  for j=0, nfrm-1 do totdn[j]=total(data1[*,*,j])
  
  sigdn=fltarr(nfrm)
  
  for j=0, nfrm-1 do sigdn[j]=robust_sigma(data1[*,*,j]/max(data1[*,*,j]))
  
  jd1=jd-jd[0]
  
  rej1=where(medimg0 gt 5000 OR totdn lt 5000, nrej1, complement=keep1) ; bad images - discard
  
  jd=jd[keep1]
  data1=data1[*,*,keep1]
  ccd_temp=ccd_temp[*,keep1]
  medimg0=medimg0[keep1]
  totdn=totdn[keep1]
  
  nfrm=n_elements(keep1)
  stop
  for j=0, nfrm-1 do begin
    
    cgimage, bytscl(data1[*,*,j], 20, 500)
    
    wait, 0.4
    
  endfor
  
  
  stop
endfor
  
  
print, 'end of program'
end