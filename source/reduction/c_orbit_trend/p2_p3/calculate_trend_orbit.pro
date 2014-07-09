pro calculate_trend_orbit

; calculate the average trend by stacking orbits according to exposure number in each sequence

  Compile_opt idl2
  !p.background=cgcolor('white')
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  ; input directory
  indir='~/BRITE/data/UB/p3/'
  
  filesin=file_search(indir+'*_p3.sav', count=nsav)
  
  obj=obj_new('IDL_Savefile', filesin[0])
  obj->restore, 'jd' 
  obj->restore, 'ccd_temp' 
  
  nfrm=n_elements(jd)
  
  roi_out=fltarr(nsav,nfrm)
  jd1=jd-jd[0]
  header_temps=ccd_temp
  roi_loc=fltarr(2,nsav)
  
  for bb=0, nsav-1 do begin
  
    obj=obj_new('IDL_Savefile', filesin[bb])
    obj->restore, 'data1'
    obj->restore, 'jd'
    obj->restore, 'xc'
    obj->restore, 'yc'
    
    roi_loc[0,bb]=xc
    roi_loc[1,bb]=yc
    
    for rr=0, nfrm-1 do begin
    
      dat=reform(data1[*,*,rr])
      
      rbin=4
      dat2=rebin_data(dat, rbin)
      
      ; find approx center of PSF
      mthr=500
      spr=2
      pks=brite_findpeaks(dat2, minthresh=mthr, spread=spr)
      
      if pks[0] eq -999.9 then goto, next_frame
      
      npks=n_elements(pks)/3.
      
      if npks gt 1 then pks=pks[*,where(pks[2,*] eq max(pks[2,*]))]
      
      pks1=round(pks)
      
      ; make a cutout centered on the target
      if pks1[0]-(8*rbin) lt 0 then x1=0 else x1=pks1[0]-(8*rbin)
      if pks1[0]+(8*rbin) gt (32*rbin)-1 then x2=(32*rbin)-1 else x2=pks1[0]+(8*rbin)
      if pks1[1]-(8*rbin) lt 0 then y1=0 else y1=pks1[1]-(8*rbin)
      if pks1[1]+(8*rbin) gt (32*rbin)-1 then y2=(32*rbin)-1 else y2=pks1[1]+(8*rbin)
      
      cutout=dat2[x1:x2,y1:y2]
      ; reorder pixels in cutout
      order1=reverse(sort(cutout))
      ncut=n_elements(order1)
      order2=order1
      order2[0]=-9999
      count=0
      snr=0
      
      ;plot, indgen(ncut), intarr(ncut), yrange=[0,100], xrange=[0,ncut], /nodata, color=cgcolor('black')
      
      repeat begin
      
        snr1=snr
        count=count+1
        
        if count eq ncut then goto, next_frame
        
        flux=total(cutout[order1[0:count]])
        
        snr=flux/sqrt(flux + count*(50)^2)
        
        order2[count]=-9999
        
        ;oplot, [count], [snr], color=cgcolor('black'), psym=2
        
      endrep until snr lt snr1
      
      if count lt 5. then goto, next_frame
      
      xx=where(order2 eq -9999)  ; xx are 1-d indices in cutout with the target flux
      ;if jd1[rr] ge 6. then begin
      ;  plot_image, cutout
      ;  oplot, (array_indices(cutout, order1[xx]))[0,*], (array_indices(cutout, order1[xx]))[1,*], color=cgcolor('orchid'), psym=2
      ;  stop
      ;endif
      
      temp_dat=dat2
      temp_dat[x1:x2,y1:y2]=-9999
      yy=where(temp_dat ne -9999, nyy)  ; xx are pixels away from target
      ;stop
      ;plot_image, dat2
      ;oplot, (array_indices(dat2, yy))[0,*], (array_indices(dat2, yy))[1,*], color=cgcolor('orchid'), psym=2
      ;print, rr
      ;wait, 1
      ; calculate the median in the in- and out- of target pixels
      roi_out[bb,rr]=robust_mean(dat2[yy],3)
      
      ;stop
      
      next_frame:
    endfor  ; end loop over individual observations
    
  endfor
  
  ; save roi_out with jd1 
    outfile='/Users/gemmawhittaker/BRITE/reduction/orbit_trends/roi_med_all/roi_med_all_p3.sav'
    save, filename=outfile, roi_out, jd1, header_temps, roi_loc
    
    print, 'End of Program'
  
end