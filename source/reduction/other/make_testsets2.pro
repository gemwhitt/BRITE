pro make_testsets2

; MODIFIED FROM make_testsets.pro on 31/3/2014 
; .... - make test sets according to alternative criteria, e.g. ccd_temp range.....

  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  
  Compile_opt idl2
  
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 0.9
  
  ; make test sets using only "good" observations - use this to characterise the data
  
  indir='~/BRITE/data/UB/p2/'
  
  outdir='~/BRITE/data/UB/testing/saturated/savfiles/'
  
  g=0
  
  cf=['*1-2','*1-7']
  
  targets=['34085','35468','39801']
  
  filesin=file_search(indir+cf[g]+'*'+targets+'*.sav', count=nsav)
  print, nsav
  
  for i=0, nsav-1 do begin
  
    fname=file_basename(filesin[i], '_p2.sav')
    
    restore, filesin[i] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, exp_time, exp_ttl, medcols
    
    
    xx=where(ccd_temp ge 16 AND ccd_temp le 18, nxx)
    
    yy=where(ccd_temp gt 24 AND ccd_temp lt 26, nyy)
    
    keep=[xx[0:59],yy[0:59]]
    
    jd=jd[keep]
    data1=data1[*,*,keep]
    ccd_temp=ccd_temp[keep]
;    exp_time=exp_time[0]
    
 
    ; get the magnitude
    vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
    
    fileout=outdir+fname+'_p2.sav'
    
    save, filename=fileout, jd, data1, ccd_temp, simbad_mags, vmag
       
  endfor
  
  
  print, 'end of program'
end

