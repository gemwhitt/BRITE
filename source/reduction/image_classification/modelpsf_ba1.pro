pro modelpsf_ba1

; PURPOSE: To devise a "model" PSF for each target based on the "best" images - get an average shape and flux - 
           ;with level of deviation for "good" images
           
; modifed from orion_modelpsf.pro on 16th June

Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.5

sat='BA'
field='ORION'

indir='~/BRITE/data/'+sat+'/p2/'+field+'/class/sav/'
  
filesin=file_search(indir+'*.sav', count=nfl)
  
plotdir='~/BRITE/data/'+sat+'/p2/'+field+'/class/models/'
  
txtout='~/BRITE/data/'+sat+'/p2/'+field+'/class/stats/psf_time_temp.txt'

rbin=20.  ; resolution for PSF-fitting
    
for ff=0, nfl-1 do begin  ; START LOOP OVER FILES
  
  restore, filesin[ff]  ; roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
                        ; medcol1, medimg0, ndead, nsat, simbad_radec, vmag, bmag, $
                        ; parlax, otype, sptype, iflag, ccd_temp2, psf_loc, flx_reg, xy_com
      
  nimg=n_elements(jd)
      
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  jd1=jd-jd[0]
  jd2=jd1[1:nimg-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits in this sequence
  gap=[-1,gap,nimg-1]
      
  ; variables to determine 
  xy_fit=fltarr(2,nimg)                     ; get x&y PSF centers for all frames
  mod1=fltarr(xdim*rbin,ydim*rbin,ngap+1)   ; get a model PSF for each orbit
  npix=intarr(ngap+1)                       ; count the number of binned PSF pixels in each orbit
  avg_temp=fltarr(ngap+1)                   ; record the average temperature in each orbit
      
  for orbit=0, ngap do begin  ; START LOOP OVER ORBITS
      
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
    
    keep=where(iflag[iloc] eq 2, nkeep, complement=rej)
    
    if nkeep lt 10 then continue  ; if there are less than 10 good frames in an orbit then skip to next
    
    iloc=iloc[keep]
        
    avg_temp[orbit]=average(ccd_temp[*,iloc])
        
    data2=data1[*,*,iloc]   ; data2 is only the images in this orbit
        
    ; DO sr-fitting to get xy_fit and models
    sr_fitting_ba1, data2,iloc,orbit,plotdir,hdname,rbin, xy_fit,mod1,npix
    
    stop
  endfor  ; END LOOP OVER ORBITS
  
  ; VIEW different model PSFs - ignore where there is NO PSF
      
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
