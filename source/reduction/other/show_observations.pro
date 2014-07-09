pro show_observations

; do plot image for every observation in sequence to verify data
; 
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

indir='~/BRITE/data/UB/p1/CENTAURUS/'

filesin=file_search(indir+'*.sav', count=nf)

for i=1, nf-1 do begin
  
  restore, filesin[i] ;jd, jd1, data1, ccd_temp, roi_dim, simbad_mags, p2_trend, med_trend
  
  nfrm=(size(data1, /dim))[2]
  
  jd1=jd-jd[0]
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  
 
  x1=670
  x2=680
  
  for j=x1, x2 DO BEGIN ;nfrm-1 do begin
    
    print, j
    
    plot_image, bytscl(data1[*,*,j], 20, 500), title=strtrim(j,2), color=cgcolor('black')
    
    wait, 1.
  endfor
  
endfor





stop
print, 'end of program'
end