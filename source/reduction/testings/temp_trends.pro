pro temp_trends

  ; Modified from map_hps
  ; Purpose: Investigate temperature trends on HPs, DC etc....
  ;
  ;
  ; CALLS: find_hps.pro
  ;
  Compile_opt idl2
  ;; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  
  hd='31237'
  
  ; input directory
  indir='~/BRITE/data/UB/p1/ORION/'
    
  filesin=file_search(indir+'*31237*.sav', count=nsav)
  fname=file_basename(filesin, '_p1.sav')
  
  if nsav eq 0 then stop
  
  fileout='~/BRITE/data/UB/test_analysis/temp_trends/plots/HD'+hd+'_temptrends.ps'
  ps_on, fileout, xsize=15, ysize=15
  
  for f=0, nsav-1 do begin
    
    hdpos=strpos(fname[0], 'HD')
    hdname=strmid(fname[0], hdpos)
  
    print, filesin[f]    
    
    ; restore data
    restore, filesin[f]; jd, jd1, data1, ccd_temp, medimg, exp_ttl

    ; number of frames in file
    nfrm=n_elements(jd)
    
    ; determine number of orbits - 1 orbit=15 mins (0.01 days), 1 gap = 85 mins (0.059 days)
    jd1=jd[1:nfrm-1]
    jdiff=jd1-jd
    gap=where(jdiff gt 0.015, ngap)  ; there are ngap+1 observations
    gap=[0,gap,nfrm-1]
    
    nhp=intarr(nfrm)  ; total number of HPs
    avhp=fltarr(nfrm) ; average value of HPs
    
    map_hps2, data1,exp_ttl, nhp,avhp
    
    ccd_temp=ccd_temp
    
    
  ;  if f eq 0 then begin
      
      plotsym, 0, /fill, 0.3
      plot, ccd_temp, nhp, psym=8, color=cgcolor('black'), xrange=[15,40], yrange=[10,120], $
        xtitle='CCD Temp (degrees C)', ytitle='Number of HPs', charsize=0.8
      ;oplot, ccd_temp, n2hp, psym=8, color=cgcolor('purple') 
      ;al_legend, ['Total num HP', 'Cluster HP'], psym=[8,8], color=[cgcolor('black'), cgcolor('purple')], /left, charsize=0.7, $
      ;  textcolor=cgcolor('black')
 ;   endif else begin 
 ;     oplot, ccd_temp, nhp, psym=8, color=cgcolor('black')
      ;oplot, ccd_temp, n2hp, psym=8, color=cgcolor('purple')
 ;   endelse
      
  ;;    if f eq 0 then plot, ccd_temp, avhp, psym=2, color=cgcolor('black'), xrange=[15,40], yrange=[0,2000] else $
  ;;    oplot, ccd_temp, avhp, psym=2, color=cgcolor('black')
    
    ; save results
    
  endfor
  
  ps_off
  
  spawn, 'open '+fileout+' &'
  stop
end
