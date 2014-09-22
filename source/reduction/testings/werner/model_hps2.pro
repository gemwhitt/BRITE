pro model_hps2

  ; modified to run over many more observations - not just one orbit

  ; program to use the PSF models to detect HPs within the PSF and correct for them
  ; Use 'modelpsf' and stacking to detect HPs in unbinned data
  ;
  Compile_opt idl2
  
  nbin=10
  
  indir='~/BRITE/TESTSETS/werner4lc/p4/'
  outdir='~/BRITE/TESTSETS/werner4lc/lc1/'
  
  filesin=file_search(indir+'*TOR*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, nf-1 do begin
  
    fname=file_basename(filesin[f], '.sav')
    
    fileout=outdir+fname+'.txt'
    
    restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1, flag, $
    ;psf_loc, npix_psf, hp_loc, nhp, modelpsf, medcol3, avg_med_diff, xy_com
    
    nimg=n_elements(jd)
    s=size(data1, /dim)
    
;    xx=where(flag ne 2, nxx)
;    for i=0, nxx-1 do begin
;      plot_image, bytscl(data1[*,*,xx[i]], 0, 200)
;      stop
;    endfor
;    stop
    jd1=jd-jd[0]
    
;    for im=0, 0 do begin
;    
;      if flag[im] ne 2 then print, 'skipped image '+strtrim(im,2)
;      if flag[im] ne 2 then continue
;      
;      mod1=modelpsf[*,*,im]
;      
;      plot_image, bytscl(mod1, 0, 500)
;      
;      plot_image, bytscl(data1[*,*,im], 0, 200)
;      
;      stop
;    endfor
    
    
    ; separate images into days
    hist1=cghistogram(jd1, binsize=1, locations=loc, reverse_indices=ri)
    ndays=n_elements(hist1)
    
    ; define new arrays to save the binned and shifted PSFs for this orbit
    binarray=fltarr(s[0]*nbin,s[1]*nbin,s[2])
    shiftarray=fltarr(s[0]*nbin,s[1]*nbin, s[2]) 
    
    for group=0, ndays-1 do begin
      
      if hist1[group] eq 0 then continue
      
      if hist1[group] lt 100 then begin
        if n_elements(iloc) eq 0 then iloc=[ri[ri[group]:ri[group+1]-1]] else $
          iloc=[iloc,ri[ri[group]:ri[group+1]-1]]
        goto, next_group
      endif
    
      if n_elements(iloc) eq 0 then iloc=ri[ri[group]:ri[group+1]-1] else $
        iloc=[iloc,ri[ri[group]:ri[group+1]-1]]
      
      ; make sure there are enough "good" images to continue with
      good=where(flag[iloc] eq 2, ngood)
      
      if ngood le 5 then begin
        flag[iloc]=0
        goto, next_group
      endif
      
      iloc=iloc[good]
      
      nim=n_elements(iloc)  ; number of images in this orbit (good ones)
      
      ; define the array for the new model..
      newmodel=fltarr(s[0]*nbin,s[1]*nbin)
      sm=size(binarray, /dim)
            
      for im=0, nim-1 do begin
        ; rebin the data using interpolation 
        im0=modelpsf[*,*,iloc[im]]  ; actual modelpsf for this image
        im1=im0*0                   ; this will be the image binned
        
        im1=congrid(im0, nbin*s[0], nbin*s[1], /interp, /center)
        
        binarray[*,*,iloc[im]]=im1
        
        ;      wset, 0
        ;      plot_image, bytscl(im0, 0, 200)
        ;
        ;      wset, 1
        ;      plot_image, bytscl(im1, 0, 200)
        ;      stop
        
        ; now modify xy_com to obtain the center of PSF in this binned image
        xcen=fix(round(xy_com[0,iloc[im]]*nbin))
        ycen=fix(round(xy_com[1,iloc[im]]*nbin))
        
        ; get the center of the binned image
        xcen2=sm[0]/2.
        ycen2=sm[1]/2.
        
        ; calculate the shifts 
        xdiff=xcen2-xcen
        ydiff=ycen2-ycen
           
        ; now shift the image to the center of the model frame
        im2=im1*0
          
        ; find the locations of ALL PSF pixels gt 0
        psfpix=where(im1 gt 0, npsf)
        xpix=(array_indices(im1,psfpix))[0,*]
        ypix=(array_indices(im1,psfpix))[1,*]
                  
        im2[xpix+xdiff,ypix+ydiff]=im1[xpix,ypix]
          
        shiftarray[*,*,im]=im2
          
        newmodel=newmodel+im2
       
      endfor
       
      newmodel=newmodel/float(nim)
      
;      plot_image, bytscl(newmodel,0, 200)
;   stop    
      ; use label-region here?
      newmod1=fltarr(sm[0]+2,sm[1]+2)
      newmod1[1,1]=newmodel
      
      thr=100
      
      ; use label region to determine number of illuminated regions and number of PSF pixels
      r2=label_region(newmod1 ge thr, /all_neighbors)
      r2=r2[1:sm[0],1:sm[1]]                         ; trim off border
      
      hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)
      
      ; sort hist2 in decreasing order
      sort2=reverse(sort(hist2))
      
      ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
      i2=array_indices(newmodel, ind2)
      
      xi=reform(i2[0,*])
      yi=reform(i2[1,*])
      
;      oplot, xi, yi, color=cgcolor('purple'), psym=2
;   stop   
      psf_aper=fltarr(sm[0],sm[1])
      psf_aper[xi,yi]=newmodel[xi,yi]
      
      ; now for this group of images get the flux inside the aperture and save to a text file
      ; with time, xcen, ycen, temp
      for im=0, nim-1 do begin
        
        im0=data1[*,*,iloc[im]]
        s=size(im0, /dim)
        
        ; bin the image
        im1=congrid(im0, nbin*s[0], nbin*s[1], /interp, /center)
        
        ; now modify xy_com to obtain the center of PSF in this binned image
        xcen=fix(round(xy_com[0,iloc[im]]*nbin))
        ycen=fix(round(xy_com[1,iloc[im]]*nbin))
        
        ; get the center of the binned image
        xcen2=s[0]*nbin/2.
        ycen2=s[1]*nbin/2.
        
        ; calculate the shifts
        xdiff=xcen2-xcen
        ydiff=ycen2-ycen
        
        ; now shift the model to the center of the PSF in the image
        mod1=psf_aper*0
        
        ; find the locations of ALL PSF pixels gt 0
        psfpix=where(psf_aper gt 0, npsf)
        xpix=(array_indices(psf_aper,psfpix))[0,*]
        ypix=(array_indices(psf_aper,psfpix))[1,*]
        
        mod1[xpix-xdiff,ypix-ydiff]=psf_aper[xpix,ypix]
        
        modpix=where(mod1 gt 0, npsf)
        xpix=(array_indices(mod1,modpix))[0,*]
        ypix=(array_indices(mod1,modpix))[1,*]
        
        flux=total(im1[xpix,ypix])
        
        temperature=average(ccd_temp[*,iloc[im]])
        
        openw, lun, fileout, /append, /get_lun
        printf, lun, jd[iloc[im]], flux, xcen/float(nbin), ycen/float(nbin), temperature, $
          format='(d14.6,x,d14.3,x,f7.3,x,f7.3,x,f7.3)'
        free_lun, lun
        
        
;        wset, 0
;        plot_image, bytscl(im1, 0, 200)
;        
;        wset, 1
;        plot_image, bytscl(mod1, 0, 200)

      
      endfor
      
;      plot_image, bytscl(psf_orbit,0, 200)
;      stop
;      for im=0, nim-1 do begin
;      
;        wset, 0
;        plot_image, bytscl(binarray[*,*,iloc[im]],0, 500)
;        
;        wset, 1
;        plot_image, bytscl(shiftarray[*,*,iloc[im]]-psf_orbit,0, 500)
;        
;        
;        stop
        undefine, iloc
        next_group:
      endfor
      
    stop
    
  endfor ; end loop over file
  
  print, 'end of program'
end