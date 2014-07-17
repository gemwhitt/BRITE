pro calculate_trend_temp

; calculate the average trend by stacking orbits according to exposure number in each sequence

  Compile_opt idl2
  !p.background=cgcolor('white')
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  ; input directory
  indir='~/BRITE/data/UB/p2/'
  
  filesin=file_search(indir+'*_p2.sav', count=nsav)
  fname=file_basename(filesin, '_p2.sav')
  
  for bb=0, nsav-1 do begin
    
    hdpos=strpos(fname[bb], 'HD')
    hdname=strmid(fname[bb], hdpos)
  
    obj=obj_new('IDL_Savefile', filesin[bb])
    obj->restore, 'data1'
    obj->restore, 'jd'
    obj->restore, 'ccd_temp'
    obj->restore, 'xc'
    obj->restore, 'yc'
    
    ; determine distance from CCD center
    ccd_xc=4008./2.
    ccd_yc=2672./2.
    
    dist1=sqrt( (ccd_xc-xc)^2. + (ccd_yc-yc)^2 )
    
    nfrm=n_elements(jd)   
    jd1=jd-jd[0]
    
    for rr=0, nfrm-1 do begin
    
      dat=reform(data1[*,*,rr])
      
      ; find approx center of PSF
      mthr=500
      spr=2
      pks=brite_findpeaks(dat, minthresh=mthr, spread=spr)
      
      if pks[0] eq -999.9 then goto, next_frame
      
      npks=n_elements(pks)/3.
      
      if npks gt 1 then pks=pks[*,where(pks[2,*] eq max(pks[2,*]))]
      
      ; make a cutout centered on the target
      if pks[0]-8 lt 0 then x1=0 else x1=pks[0]-8
      if pks[0]+8 gt 31 then x2=31 else x2=pks[0]+8
      if pks[1]-8 lt 0 then y1=0 else y1=pks[1]-8
      if pks[1]+8 gt 31 then y2=31 else y2=pks[1]+8
      
      cutout=dat[x1:x2,y1:y2]
      ; reorder pixels in cutout
      order1=reverse(sort(cutout))
      ncut=n_elements(order1)
      order2=order1
      order2[0]=-9999
      count=0
      snr=0
      
      repeat begin
      
        snr1=snr
        count=count+1
        
        if count eq ncut then goto, next_frame
        
        flux=total(cutout[order1[0:count]])
        
        snr=flux/sqrt(flux + count*(30)^2)
        
        order2[count]=-9999
        
      endrep until snr lt snr1
      
      if count lt 5. OR count gt float(ncut)/2. then goto, next_frame
      
      xx=where(order2 eq -9999)  ; xx are 1-d indices in cutout with the target flux
      if jd1[rr] ge 6. then begin
        plot_image, cutout
        oplot, (array_indices(cutout, order1[xx]))[0,*], (array_indices(cutout, order1[xx]))[1,*], color=cgcolor('orchid'), psym=2
        stop
      endif
      
      temp_dat=dat
      temp_dat[x1:x2,y1:y2]=-9999
      yy=where(temp_dat ne -9999, nyy)  ; xx are pixels away from target
      ;stop
      ;plot_image, dat
      ;oplot, (array_indices(dat, yy))[0,*], (array_indices(dat, yy))[1,*], color=cgcolor('orchid'), psym=2
      ;stop
      
      ; caulculate the median in the in- and out- of target pixels
      if n_elements(in_pix) eq 0 then begin
        in_pix=median(cutout[order1[xx]])
        out_pix=median(dat[yy])
        temp1=ccd_temp[rr]
      endif else begin
        in_pix=[in_pix,median(cutout[order1[xx]])]
        out_pix=[out_pix, median(dat[yy])]
        temp1=[temp1,ccd_temp[rr]]
      endelse
      
      next_frame:
    endfor  ; end loop over individual observations
    
    result=histogram(temp1, binsize=0.5, locations=loc, reverse_indices=ri)
    
    ;window, 0, xpos=1500, ypos=200
    ;cghistoplot, temp1, binsize=0.5, backcolorname='white', datacolorname='blue'
    
    
    med_dn=fltarr(n_elements(result))
    for ii=0, n_elements(result)-1 do med_dn[ii]=median(out_pix[ri[ri[ii]:ri[ii+1]-1]])
    
    rob_dn=fltarr(n_elements(result))
    for ii=0, n_elements(result)-1 do mean_dn[ii]=robust_mean(out_pix[ri[ri[ii]:ri[ii+1]-1]],2)
    
    ;window, 1, xpos=2300, ypos=200
    if bb eq 0 then plot, loc, mean_dn, xrange=[15,28], xstyle=1, yrange=[80,110], ystyle=1, color=cgcolor('black'), psym=2, /nodata
    if dist1 gt 700 then oplot, loc, mean_dn, color=cgcolor('orchid'),psym=2
    if dist1 lt 600 then oplot, loc, mean_dn, color=cgcolor('green'),psym=2
    
    
    print, hdname
    print, robust_mean(rob_dn,2)
    stop
    
    
    
    
  endfor
  
  stop
  
end