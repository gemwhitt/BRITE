pro psf_aper

  ; Use modelpsf AND psf_loc to build an optimised aperture to apply to all images to get photometry
  ; This is in binned image mode
  ; 
  ; Output: Save p5 files with the binned data and aperture info
  ;
  Compile_opt idl2
  
  nbin=10
  
  sat='BA'
  
  field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p4/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/p5/'
  
  chkout=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
  
  filesin=file_search(indir+'*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, 2 do begin ;nf-1 do begin
  
    fname=file_basename(filesin[f], '.sav')
        
    restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag, nhp1,  $
      ;psf_loc, npix_psf, modelpsf, xy_com
    
    nimg=n_elements(jd)
    s=size(data1, /dim)
    
    jd1=jd-jd[0]
    
    good=where(flag eq 2, ngood)
    
    ; define new arrays to save the binned and shifted PSFs for this orbit
    binarray=fltarr(s[0]*nbin,s[1]*nbin,ngood)
    shiftarray=fltarr(s[0]*nbin,s[1]*nbin,ngood)
    pos=-1
        
    ; define the array for the averaged psf
    newmodel=fltarr(s[0]*nbin,s[1]*nbin)
    sb=size(newmodel, /dim)
      
    for im=0, nimg-1 do begin                ; loop over each image, shift the modelpsf to the center and stack the data
      
      if flag[im] ne 2 then continue
      
      pos=pos+1
      
      ; rebin the data using interpolation
      im0=modelpsf[*,*,im]        ; actual modelpsf for this image
        
      im1=congrid(im0, nbin*s[0], nbin*s[1], /interp, /center)  ; the binned image
        
      binarray[*,*,pos]=im1
        
;            wset, 0
;            plot_image, bytscl(im0, 0, 200)
;      
;            wset, 1
;            plot_image, bytscl(im1, 0, 200)
;            stop
        
        ; now modify xy_com to obtain the center of PSF in this binned image
        xcen=fix(round(xy_com[0,im]*nbin))
        ycen=fix(round(xy_com[1,im]*nbin))
        
        ; get the center of the binned image
        xcen2=sb[0]/2.
        ycen2=sb[1]/2.
        
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
        
;        wset, 0
;        plot_image, bytscl(im1, 0, 200)
;        
;        wset, 1
;        plot_image, bytscl(im2, 0, 200)
;        stop
;        
        shiftarray[*,*,pos]=im2
        
        newmodel=newmodel+im2
        
      endfor
      
      newmodel=newmodel/float(nimg)
      
      plot_image, bytscl(newmodel,0, 500)
      
      sos_aper, newmodel, nbin, sos_aper1
      
      ; Old method.....      
;      thr=200
;      
;      ; use label region to determine number of illuminated regions and number of PSF pixels
;      r2=label_region(newmod1 ge thr, /all_neighbors)
;      r2=r2[1:sb[0],1:sb[1]]                         ; trim off border
;      
;      hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)
;      
;      ; sort hist2 in decreasing order
;      sort2=reverse(sort(hist2))
;      
;      ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
;      i2=array_indices(newmodel, ind2)
;      
;      xi=reform(i2[0,*])
;      yi=reform(i2[1,*])
;      
;;      oplot, xi, yi, color=cgcolor('purple'), psym=2
;;      stop
;       
;      psf_aper=fltarr(sb[0],sb[1])
;      psf_aper[xi,yi]=newmodel[xi,yi]

psf_aper=sos_aper1
      
      ; save p5 file with binned data and psf_aper
      fileout=outdir+file_basename(filesin[f])
      save, filename=fileout, roi_name, jd, data1, roi_loc, ccd_temp, $
        vmag, bmag, flag, psf_loc, npix_psf, modelpsf, xy_com, psf_aper
      
  endfor ; end loop over file
  
  print, 'end of program'
  print, 'Now do aper_phot'
  
end