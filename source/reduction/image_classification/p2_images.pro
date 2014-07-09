pro p2_images

  ; check p2 images try to establish trends and thresholds for image categories...
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 1.5
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  ; use test data
  ;indir='~/BRITE/data/UB/testing/ORION/p1_subset1/'
  indir='~/BRITE/data/UB/p2/ORION/'
  
  outdir='~/BRITE/data/UB/plots/ORION/systematic_trends/'
  
  filesin=file_search(indir+'*.sav', count=nf)
  
  for i=0, nf-1 do begin
  
    restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
                        ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype, iflag
    
    good=where(iflag eq 2, ngood, complement=bad)
    
    nfrm=ngood
    
    xdim=(size(data1, /dim))[0]
    ydim=(size(data1, /dim))[1]
    
    data1=data1[*,*,good]
    jd=jd[good]
    ccd_temp=ccd_temp[*,good]
    ccd_temp2=fltarr(nfrm)
    
    for jj=0, nfrm-1 do ccd_temp2[jj]=avg(ccd_temp[*,jj])
    
    ccd_temp=ccd_temp2
    
    jd1=jd-jd[0]
    
    npix=[]
    flx=[]
    
    for j=0, nfrm-1 do begin
      
      data2=reform(data1[*,*,j])
      
      totdn=total(data2)
      
      mostdn=0.5*totdn
      
      data3=reform(data2, xdim*ydim)
      sort1=reverse(sort(data3))
      data3=data3[sort1]
      subtot=fltarr(xdim*ydim)
      for jj=0, xdim*ydim-1 do subtot[jj]=total(data3[0:jj])
      xx=(where(subtot ge mostdn))[0]
      
      ; save xx = number of pixels containing most of the flux + save what most flux is!
      npix=[npix,xx]
      flx=[flx,mostdn]
              
    endfor
            
    sub1=where(npix le 50, nsub1)
    med1=median(flx[sub1])
    sig1=robust_sigma(flx[sub1])
    
    ;cghistoplot, flx[sub1], binsize=500, locations=loc, histdata=result
    ;oplot, [med1-2*sig1,med1-2*sig1], [0, 20000], color=cgcolor('blue') 
    ;stop
    
    keep1=where(flx[sub1] ge med1-2*sig1 AND flx[sub1] lt med1+2*sig1, nkeep1, complement=rej1)
    
    sg1=sub1[keep1]
    nsg1=nkeep1
    
    xx=where(jd1 ge 151)
    
    for ii=0, 1000 do begin
    
     plot_image, bytscl(data1[*,*,35115+ii], 20, 100), title=ccd_temp[35115+ii], color=cgcolor('black')
     wait, 0.3
     
   endfor
   
   stop
    
    plot, jd1, ccd_temp, color=cgcolor('black'), psym=2, xtitle='Time (days)', ytitle='Temperature', $
      title=roi_name[0], charsize=1.
    oplot, jd1[sg1], ccd_temp[sg1], color=cgcolor('purple'), psym=2
    
    
    stop
  endfor
  
  
  print, 'end of program'
end