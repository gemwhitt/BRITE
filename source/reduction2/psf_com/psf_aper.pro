pro psf_aper, sat, field, target

  ; Use modelpsf AND psf_loc to build an optimised aperture to apply to all images to get photometry
  ; This is in binned image mode
  ; 
  ; Output: Save p5 files with the binned data and aperture info
  ;
  Compile_opt idl2
  
  nbin=10
  
;  sat='BA'
  
;  field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p4/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/p5/'
  
  chkout=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
    
  filesin=file_search(indir+target+'*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for f=18, nf-1 do begin
  
    fname=file_basename(filesin[f], '_p4.sav')
        
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
      
      ;sos_aper, newmodel, nbin, sos_aper1
      
      thr_aper, newmodel, nbin, thr_aper1
      
    
psf_aper=thr_aper1

;plot_image, bytscl(psf_aper[*,*,1],0, 500)
;stop
      
      ; save p5 file with binned data and psf_aper
      fileout=outdir+fname+'.sav'
      
      save, filename=fileout, roi_name, jd, data1, roi_loc, ccd_temp, $
        vmag, bmag, flag, psf_loc, npix_psf, modelpsf, xy_com, psf_aper, medimg0, medcol
      
  endfor ; end loop over file
  
  print, 'end of program'
  print, 'Now do aper_phot'
  
end