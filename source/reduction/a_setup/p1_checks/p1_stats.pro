pro p1_stats

; created 11/04/2014 - similar to p1_info - test on Centaurus data
;
; get date range, plot dates, ccd temps, median image ....
;
Compile_opt idl2

!p.background=cgcolor('white')
plotsym, 0, 0.9, /fill

nstk=1
  
;indir='~/BRITE/data/UB/p1/ORION/'
indir='~/BRITE/data/UB/p2/CENTAURUS/'+strtrim(nstk,2)+'stk/'

filesin=file_search(indir+'*.sav', count=nsav)

outdir='~/BRITE/data/UB/plots/CENTAURUS/systematic_trends/medimg_2_'+strtrim(nstk,2)+'/'
  
for i=0, nsav-1 do begin
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                      ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
                      
  hdname=roi_name[0]                   
                      
  nfrm=n_elements(jd)
            
  jd1=jd-jd[0] 
  
  medimg2=fltarr(nfrm)
  for j=0, nfrm-1 do medimg2[j]=median(data1[*,*,j])
  
  fileout=outdir+hdname+'_'+strtrim(nstk,2)+'.ps'
  
;  ps_on, fileout, xsize=18, ysize=28

;  !p.multi=[0,1,2,0,0]               
  cgplot, ccd_temp, medimg, color='black', psym=8, xtitle='CCD Temp (celsius)', ytitle='Image Median (DN)', $
    charsize=0.7, xstyle=1, xrange=[28,42], title=hdname, yrange=[0,20000]
  
;  cgplot, ccd_temp, medimg2, color='green', psym=8, charsize=0.7, xstyle=1, xrange=[28,42], title='P2 Image Median', /overplot
    
  ; get number of orbits
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  gap=[0,gap,nfrm-1]

  ; loop over orbits and plot
  for j=0, ngap do begin
    
    if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+(gap[j]+1)
    
    if j eq 0 then cgplot, ccd_temp[iloc], medimg[iloc], color='black', psym=8, xtitle='Points in orbit', $
      ytitle='Image Median', charsize=0.7, yrange=[80,300], xrange=[30,40] else $
      cgplot, ccd_temp[iloc], medimg[iloc], color='black', psym=8, /overplot
      
   ; cgplot, ccd_temp[iloc], medimg[iloc], color='black', psym=8, xtitle='Points in orbit', $
   ;   ytitle='Image Median', charsize=0.9, xrange=[30,45]
      
      stop
   
  endfor

 ; ps_off
 
endfor
  
print, 'end of program'
end