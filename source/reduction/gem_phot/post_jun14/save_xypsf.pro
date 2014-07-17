pro save_xypsf

  ; Main program to run get_xypsf.pro and save results in each .sav file
  ; to be called by p4_gem7.pro
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, 0.8, /fill
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  indir='~/BRITE/data/UB/p2/ORION/class/sav/'
  
  filesin=file_search(indir+'*.sav', count=nfiles)
  
  fname=file_basename(filesin, '_p2.sav')
  
  outdir=indir
  
  for ff=0, 0 do begin  ;nfiles-1 do begin
  
    print, fname[ff]
    
    restore, filesin[ff]  ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
                          ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
                          ;parlax, otype, sptype, iflag, ccd_temp2
    
    hdname=fname[ff]
    
    
    ;get total number of frames in this file
    nfrm=n_elements(jd)
    
    jd1=jd-jd[0]
    jd2=jd1[1:nfrm-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap)
    gap=[-1,gap,nfrm-1]
    
    max_dn=lonarr(nfrm)
    
    for i=0, nfrm-1 do max_dn[i]=max(data1[*,*,i])
    
    xy_psf=fltarr(2,nfrm)
    
    for g=0, ngap do begin
  
      iloc=indgen(gap[g+1]-gap[g])+gap[g]+1
      
      img_cl=iflag[iloc]
      
      good_im=where(img_cl eq 2, ngood, complement=bad_im)
      
      if ngood lt 2 then continue else iloc=iloc[good_im]
    
      sdat=data1[*,*,iloc]  ; subset of dat - for one orbit
    
      ; get the PSF centers
      sr_fit_xypsf, sdat,iloc,g, xy_psf
    
    endfor
    
  ;  good_pos=where(xy_psf[0,*] ne 0.0 AND xy_psf[1,*] ne 0.0, ngood, complement=bad_pos)
    
    ;plot, xy_psf[0,good_pos], max_dn[good_pos], color=cgcolor('black'), psym=2

  ; re-save .sav file with xy_psf
  fileout=filesin[ff]  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
  medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
  parlax, otype, sptype, iflag, ccd_temp2, xy_psf
    
  endfor
  
  print, 'end of program'
end

