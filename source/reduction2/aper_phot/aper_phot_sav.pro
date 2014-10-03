pro aper_phot_sav, sat, field, target

  ; make a sav file with the results of aper_phot rather than a txt file 

  ; Extract photometry from images using aper_psf and xy_com with data1
  
  ; Save to a text file with all variables
  
  ; 2 apertures in psf_aper
  
  Compile_opt idl2
  
  nbin=10
  
 ; sat='BA'
  
 ; field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p5/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_sav/'
  
  chkout=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
    
  filesin=file_search(indir+target+'*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, nf-1 do begin
  
    fname=file_basename(filesin[f], '.sav')
    
    obj=obj_new('IDL_Savefile', filesin[f])
    obj->restore, 'jd'
    obj->restore, 'data1'
    obj->restore, 'ccd_temp'
    obj->restore, 'vmag'
    obj->restore, 'bmag'
    obj->restore, 'flag'
    obj->restore, 'psf_loc'
    obj->restore, 'xy_com'
    obj->restore, 'psf_aper'
    obj->restore, 'medimg0'
    obj->restore, 'medcol'
    destroy, obj
    
    nimg=n_elements(jd)
    
    nap=(size(psf_aper, /dim))[2]
    
    s=size(data1, /dim)
    
    jd1=jd-jd[0]
    
    good=where(flag eq 2, ngood)
    
    ; if there is already an output file for this target, then delete it
    chk=file_search(outdir+fname+'.txt', count=nchk)
    if nchk gt 0 then spawn, 'rm '+outdir+fname+'.txt'
    
    flux=fltarr(ngood,nap)
    resid=fltarr(ngood,nap)
    time=[]
    xcom=[]
    ycom=[]
    temperature=[]
    frame=good
    bkgd=[]
    bkgd_err=[]
    
    
    for im=0, ngood-1 do begin
    
      ind=frame[im]
      
      ; define new arrays to save the binned and shifted PSFs for this orbit
      binarray=fltarr(s[0]*nbin,s[1]*nbin)
      shiftarray=fltarr(s[0]*nbin,s[1]*nbin)
      
      im0=data1[*,*,ind]
      
      ; calulcate the background using psf_loc (to avoid the PSF pixels)
      backimg=float(im0)
      backimg[psf_loc[0,ind]:psf_loc[1,ind],psf_loc[2,ind]:psf_loc[3,ind]]=!Values.F_NAN
      
      backgd=backimg[where(finite(backimg) eq 1, nbkgd)]
      
      fr = fractile(backgd, [.25,.75])
      
      xx=where(backgd ge fr[0] AND backgd le fr[1], nxx)
      
      ;   print, total(backgd[xx])/float(nxx)
      
      bkgd=[bkgd,total(backgd[xx])/float(nxx)]
      
      bkgd_err=[bkgd_err,sqrt((total((backgd[xx]-bkgd)^2))/float(nxx))]
      
      
      
      ; HISTOGRAM METHOD.....
      ; do histogram of background pixels and pick the value which is most represented but not zero
      ;    result=histogram(backgd, min=1, binsize=1, locations=loc)
      ;
      ;    sort1=reverse(sort(result))
      ;
      ;    ; choose the most populated value
      ;    bkgd=loc[sort1[0]]
      ;    bkgd_err=sqrt((total((backgd-mean(backgd))^2))/float(nbkgd))
      
      ; AVERAGE method
      
      
      ; remove the background value from im0
      im0=im0-bkgd[im]
      
      ; bin the image
      im1=congrid(im0, nbin*s[0], nbin*s[1], /interp, /center)
      
      ; now modify xy_com to obtain the center of PSF in this binned image
      xcen=fix(round(xy_com[0,ind]*nbin))
      ycen=fix(round(xy_com[1,ind]*nbin))
      
      ; get the center of the binned image
      xcen2=s[0]*nbin/2.
      ycen2=s[1]*nbin/2.
      
      ; calculate the shifts
      xdiff=xcen2-xcen
      ydiff=ycen2-ycen
      
      ; now shift aperture to the center of the PSF in the image
      sap=size(psf_aper, /dim)
      nap=sap[2]
      
      for ap=0, nap-1 do begin
        modelapp=reform(psf_aper[*,*,ap])
        
        aper=modelapp*0
        
        ; find the locations of ALL PSF pixels gt 0
        psfpix=where(modelapp gt 0, npsf)
        xpix=(array_indices(modelapp,psfpix))[0,*]
        ypix=(array_indices(modelapp,psfpix))[1,*]
        
        aper[xpix-xdiff,ypix-ydiff]=modelapp[xpix,ypix]
        
        modpix=where(aper gt 0, npsf)
        xpix=(array_indices(aper,modpix))[0,*]
        ypix=(array_indices(aper,modpix))[1,*]
        
        flux[im,ap]=total(im1[xpix,ypix])/nbin^2.
        
;        plot_image, bytscl(im1, 0, 500)
;        oplot, xpix, ypix, color=cgcolor('purple'), psym=2
;        stop
      
        ; calculate the residual in the image after subtracting the flux within the model aperture
        im2=im1
        im2[xpix,ypix]=0
        resid[im,ap]=total(im2)/nbin^2.
        
      endfor
      
      time=[time,jd[ind]]
      xcom=[xcom,xy_com[0,ind]]
      ycom=[ycom,xy_com[1,ind]]
      temperature=[temperature,average(ccd_temp[*,ind])]
      
    endfor  ; end loop over image
    
    medimg=medimg0[frame]
    data=data1[*,*,frame]
    
    
;    plotsym, 0, /fill, 0.2
;    plot, time, flux[*,0], color=cgcolor('black'), psym=8
;    stop
;    
    fileout=outdir+fname+'.sav'
    save, filename=fileout, time, frame, flux, bkgd, bkgd_err, xcom, ycom, temperature, $
      vmag, bmag, resid, medimg, psf_aper, data
 
    
    
  endfor  ; end loop over file
  
  print, 'End of Program'
  print, 'Now do sigma clipping with resid, then interpix, then sigma clip flux using smooth'
  
end

