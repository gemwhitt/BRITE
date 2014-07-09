pro test_sensitivity
  
  ;sunday january 17th
  ;
  Compile_opt idl2
  
  ;
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  plotsym, 0, 0.8, /fill
  ;
  ;
  indir='~/BRITE/data/UB/p3_2/'
  filesin=file_search(indir+'Orion*.sav', count=nfiles)
  fname=file_basename(filesin, '_p3.sav')
  
  outdir='~/BRITE/data/UB/sensitivity/'
  
  for bb=0, nfiles-1 do begin
    ;stop
    print, fname[bb]
    
    restore, filesin[bb] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
    ; print the magnitude
    vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
    print, 'vmag='+vmag
    
    dat=data1 ; no effect to data1
    
    ; get total number of frames in this file
    nfrm=(size(dat, /dim))[2]
    
    max_dn=fltarr(nfrm)
      
    for img=0, nfrm-1 do begin
       
      dat1=dat[*,*,img]
                  
      max_dn[img]=max(dat1)
      
      if img eq 0 OR img eq (nfrm-1) then begin
        plot_image, bytscl(dat1, 50, 500)
        stop
      endif
      
    endfor ; END LOOP OVER THIS SUB-SECTION OF ORBITS - IS USING
    
    jd=jd-jd[0]
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  new code 15th december
    xx=where(jd-jd[0] ge 0.0 AND jd lt 1.8, nxx)
    
    max_dn=max_dn[xx]
    ccd_temp=ccd_temp[xx]
    
    stop 
  endfor  ; end loop over file
  
 
  print, 'end of program'
  
end