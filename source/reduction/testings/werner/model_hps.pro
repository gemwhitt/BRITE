pro model_hps

; program to use the PSF models to detect HPs within the PSF and correct for them
; Use 'modelpsf' and stacking to detect HPs in unbinned data
; 
Compile_opt idl2

nbin=10

indir='~/BRITE/TESTSETS/werner4lc/p4/'

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin

  fname=file_basename(filesin[f], '.sav')
  
  restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
  ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1, flag, $
  ;psf_loc, npix_psf, hp_loc, nhp, modelpsf, medcol3, avg_med_diff, xy_com
  
  nimg=n_elements(jd)
  s=size(data1, /dim)
  
  ; define new arrays to save the binned and shifted PSFs for this orbit
  binarray=fltarr(s[0]*nbin,s[1]*nbin,s[2])
  shiftarray=fltarr(s[0]*nbin,s[1]*nbin, s[2])
  
  jd1=jd-jd[0]
  
    for im=0, 0 do begin
  
      if flag[im] ne 2 then print, 'skipped image '+strtrim(im,2)
      if flag[im] ne 2 then continue
  
      mod1=modelpsf[*,*,im]
  
      plot_image, bytscl(mod1, 0, 500)
  
      wait, 0.5
    endfor
  stop  
  
  ; separate images into days
  jd2=jd1[1:nimg-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  gap=[-1,gap,nimg-1]
  
  for orbit=0, ngap do begin
    
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
    
    ; make sure there are enough "good" images to continue with
    good=where(flag[iloc] eq 2, ngood)
    
    if ngood le 5 then begin
      flag[iloc]=0
      goto, next_orbit
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
        xcen2=sm[0]*nbin/2.
        ycen2=sm[1]*nbin/2.
        
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
          
        newmodel=newmodel+im1
       
      endfor
       
    
    newmodel=newmodel/float(nim)
    
    plot_image, bytscl(newmodel,0, 200)
  stop  
    sm=size(newmodel, /dim)
    
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
    
    oplot, xi, yi, color=cgcolor('purple'), psym=2
    
    psf_orbit=fltarr(sm[0],sm[1])
    psf_orbit[xi,yi]=newmodel[xi,yi]
    
    plot_image, bytscl(psf_orbit,0, 200)
    stop
    for im=0, nim-1 do begin
      
    wset, 0
    plot_image, bytscl(data2[*,*,im],0, 500)
    
    wset, 1
    plot_image, bytscl(data2[*,*,im]-psf_orbit,0, 500)
    
      
    stop
    
    endfor
      
    
    
   
   next_orbit: 
  endfor
  stop

endfor ; end loop over file

print, 'end of program'
end