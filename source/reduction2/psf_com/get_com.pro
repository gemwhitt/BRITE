pro get_com, sat, field, target

  ; Purpose: Calulate the COM for each target PSF and save info in the .sav file
  
  Compile_opt idl2
  
 ; sat='BA'
  
 ; field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p4/'
  
  filesin=file_search(indir+target+'*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, nf-1 do begin
  
    fname=file_basename(filesin[f], '.sav')
    
    restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag, nhp1,  $
    ;psf_loc, npix_psf, modelpsf
    
    nimg=n_elements(jd)
    
    jd1=jd-jd[0]
    
    xy_com=fltarr(2,nimg)
    
    for im=0, nimg-1 do begin ; begin loop over images
    
      if flag[im] ne 2 then continue
      
      im0=reform(data1[*,*,im])
      
      s0=size(im0, /dim)
      
      ; Use psf_loc to get cutout
      x1=(psf_loc[0,im]-1 > 0)
      x2=(psf_loc[1,im]+1 < s0[0]-1)
      y1=(psf_loc[2,im]-1 > 0)
      y2=(psf_loc[3,im]+1 < s0[1]-1)
      
      if (x2-x1)*(y2-y1) lt 25 OR (x2-x1)*(y2-y1) gt 0.33*s0[0]*s0[1] then begin
        flag[im]=0
        goto, next_im
      endif
      
      cutout=im0[x1:x2,y1:y2]
      
      s=size(cutout, /dim) ; length of array in x and y
      
      imtot=total(cutout)
      
      x_cen=total(total(cutout, 2)*indgen(s[0]))/imtot
      y_cen=total(total(cutout, 1)*indgen(s[1]))/imtot
      
      xy_com[0,im]=x_cen+x1
      xy_com[1,im]=y_cen+y1
      
       ;plot results to check
;            wset, 0
;            plot_image, bytscl(modelpsf[*,*,im], 20, 200)
;            oplot, [xy_com[0,im]], [xy_com[1,im]], color=cgcolor('purple'), psym=2
;            stop
      
      next_im:
    endfor  ; end loop over images
    
    
    ; re-save .save file with xy_com added
    
    save, filename=filesin[f], roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag,  $
      psf_loc, npix_psf, modelpsf, xy_com
      
  endfor ; end loop over file
  
  print, 'end of program'
end