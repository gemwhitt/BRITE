pro orion_modelpsf

  ; modified from orion_psf_sg1
  ; use sr_fitting_sg1.pro
  ; save xy_psf centers for each image
  ; save a model PSF and average flux
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 1.5
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  indir='~/BRITE/data/UB/p2/ORION/class/sav/'
  
  filesin=file_search(indir+'*.sav', count=nfl)
  
  outdir='~/BRITE/data/UB/p2/ORION/class/plots/'
  
  txtout='~/BRITE/data/UB/p2/ORION/class/stats/psf_time_temp.txt'
  
  sdate=[0,137,149.5]
  stemp=[18,28,18]
  
  for ff=6, nfl-1 do begin
  
    restore, filesin[ff]  ; roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
                          ; medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
                          ; parlax, otype, sptype, iflag, ccd_temp2
                        
    xypsf2=xy_psf
    
    ; select files < 10 days
    jd1=jd-jd[0]
    
    fluxes=dblarr(3)
    npixels=fltarr(3)
    
    for kk=1, 1 do begin
         
    xx=where(jd1 ge sdate[kk] AND jd1 le sdate[kk]+10, nxx)
    
    jd0=jd1[xx]
    data0=data1[*,*,xx]
    ccd_temp0=ccd_temp2[xx]
    iflag0=iflag[xx]
    
    nfrm=nxx
       
    xdim=(size(data0, /dim))[0]
    ydim=(size(data0, /dim))[1]
    
    jd2=jd0[1:nfrm-1]
    jdiff=jd2-jd0
    gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits in this sequence
    gap=[-1,gap,nfrm-1]
    
    ; variables to determine in sr_fitting_sg1
    rbin=10.
    xy_psf=fltarr(2,nfrm) ; for all frames
    mod1=fltarr(xdim*rbin,ydim*rbin,ngap+1)
    fl=dblarr(ngap+1)
    npix=intarr(ngap+1)
    avg_temp=fltarr(ngap+1)
        
    for orbit=0, ngap do begin
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      
      class=iflag0[iloc]
      bad_img=where(class ne 2, nbad, complement=good_img)
      
      if n_elements(good_img) lt 20 then continue
            
      avg_temp[orbit]=median(ccd_temp0[iloc])
      
      data2=data0[*,*,iloc]
      
      nimg=n_elements(iloc)
      
      plot_image, bytscl(data2[*,*,20], 50, 500)
      stop
     
      ; get the PSF centers
      sr_fitting_sg1, data2,iloc,orbit,outdir,hdname, xy_psf,mod1,npix
      
    endfor
    
    xx=where(avg_temp ge stemp[kk] AND avg_temp le stemp[kk]+4, nxx)
    
    model=reform(mod1[*,*,0])*0.0
    
    for jj=0, nxx-1 do model=model+mod1[*,*,xx[jj]]
    
    model0=model/float(nxx)  
    modpix=where(model0 gt 0, npix)
    
    fluxes[kk]=total(model0)/(rbin^2)
    npixels[kk]=npix
   
    ; plot image
    ;ps_on, outdir+roi_name[0]+'_model'+strtrim(kk,2)+'.ps', xsize=15, ysize=14
    plot_image, bytscl(model0, 50, 500);, $
    ;  title=roi_name[0]+', '+strtrim(sdate[kk],2)+' days, '+strtrim(stemp[kk]+2.,2)+' deg', $
    ;  charsize=0.8, color=cgcolor('black')
    ;  ps_off
    stop
    endfor
    
    ; save output
   ; openw, lun, txtout, /get_lun, /append
   ; printf, lun, roi_name[0], vmag, fluxes, npixels, $
   ;   format='(a10,x,f7.2,x,d14.1,x,d14.1,x,d14.1,x,f7.1,x,f7.1,x,f7.1)'
   ; free_lun, lun
   
   xy_psf=xypsf2
   ; save the model in the p2 file
   save, filename=filesin[ff], roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
    medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
    parlax, otype, sptype, iflag, ccd_temp2, model0, xy_psf
   
   
   
   continue
    
    ; get frames with iflag eq 2 within jd1 < 10.
    tar_frm=where(iflag eq 2, ntar)    
    ratio_flx=fltarr(ntar)
    
    for im=0, ntar-1 do begin
    
      data2=reform(data1[*,*,tar_frm[im]])
      data2=rebin_data(data2, rbin)
     
      xc=xy_psf[0,tar_frm[im]]
      yc=xy_psf[1,tar_frm[im]]
      
      ; shift the model to match the center of the image PSF
      ccd_xc=(size(data2, /dim))[0]/2.
      ccd_yc=(size(data2, /dim))[1]/2.
      
      mod2=model*0.0 ; new model - shifted
      ; get index locations of model PSF in model2 - ixpsf and iypsf
      ixpsf=(array_indices(model, where(model gt 0)))[0,*]
      iypsf=(array_indices(model, where(model gt 0)))[1,*]
      shx=xc-ccd_xc
      shy=yc-ccd_yc
      
      mod2[ixpsf+shx,iypsf+shy]=model[ixpsf,iypsf]
      
      ; get total flux in image within PSF pixels
      im_flx=total(data2[ixpsf+shx,iypsf+shy])/(rbin^2)
           
      res=data2-mod2
      
      ratio_flx[im]=im_flx/mod_flx*100.
          
    endfor
    
  endfor
  
  print, 'end of program'
end
