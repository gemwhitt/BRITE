pro get_psf_boundary, sat, field, target

  ; Program to determine frames with an adequate number of PSF pixels - use label_region.pro
  ; 
  ; Use p4 images - which have cleaned of HPs
  ;
  ; Output: Save a new variation of the input files, but in the date range of TOR and with flag added to denote
  ; good/bad/blank frames
  ;
  Compile_opt idl2
  
  ;sat='BA'
  
  ;field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p4/'
  
  filein=file_search(indir+target+'*.sav', count=nf)
  
  outdir=indir
  
  for f=0, nf-1 do begin
  
    print, file_basename(filein[f], '.sav')
    
    restore, filein[f]  ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag, nhp1
    
    nimg=n_elements(jd)
    
    ; establish psf_loc array for recording the boundary of the PSF pixels - use to create a model?
    psf_loc=roi_loc*0
    npix_psf=intarr(nimg)
    ; record exact PSF pixels and their values - save in separate array
    modelpsf=data1*0.
    
    jd1=jd-jd[0]
    
    ; set threshold for finding illuminated pixels
    thr=50
    
    for im=0, nimg-1 do begin
    
      if flag[im] eq 0 then goto, skipimage
      
      dat=data1[*,*,im]
      
      s=size(dat, /dim)
      
      ; add a border to the image
      dat2=lonarr(s[0]+2,s[1]+2)
      dat2[1,1]=dat
      
      ; use label region to determine number of illuminated regions and number of PSF pixels
      r2=label_region(dat2 ge thr, /all_neighbors)
      r2=r2[1:s[0],1:s[1]]                         ; trim off border
      
      hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)
      
      ; sort hist2 in decreasing order
      sort2=reverse(sort(hist2))
      
      if n_elements(hist2) eq 1 then begin  ; no target OR bad
        flag[im]=0
        goto, skipimage
      endif
      
       ;the number of results in hist2 are the number of regions (blobs) > thr, plus the background pixels
       ;so if n_elements(hist2) > 2 then - HP!
      if n_elements(hist2) gt 2 then begin
                
        ; record locations of HPs
        for region=2, n_elements(hist2)-1 do begin
          ; find 1d locations of HPs in each section
          ihp=ri2[ri2[sort2[region]]:ri2[sort2[region]+1]-1]
          i2hp=array_indices(dat, ihp)
          
          for hp=0, n_elements(ihp)-1 do begin
            
            x1=(i2hp[0,hp]-1) > 0
            x2=(i2hp[0,hp]+1) < s[0]-1
            y1=(i2hp[1,hp]-1) > 0
            y2=(i2hp[1,hp]+1) < s[1]-1
            
           
            dat[i2hp[0,hp], i2hp[1,hp]]=median(dat[x1:x2,y1:y2])
          
          
           data1[*,*,im]=dat
          endfor
         
          ;          plot_image, bytscl(dat, 0, 200)
          ;          oplot, i2hp[0,*], i2hp[1,*], color=cgcolor('orange'), psym=2         
        endfor
      endif
      
      ; now get 2nd biggest group of pixels i.e. PSF pixels
      ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
      i2=array_indices(dat, ind2)
      
      xi=reform(i2[0,*])
      yi=reform(i2[1,*])
      
      ; get borders of PSF
      psf_loc[*,im]=[min(xi),max(xi),min(yi),max(yi)]
      
      ;      wset, 0
      ;      plot_image, bytscl(dat, 0, 200)
      ;      ;oplot, i2[0,*], i2[1,*], color=cgcolor('green'), psym=2
      
      npix_psf[im]=hist2[sort2[1]]
      
      ; record PSF for model
      mod1=dat*0
      mod1[xi,yi]=dat[xi,yi]
      
      modelpsf[*,*,im]=mod1
      
      ;      wset, 1
      ;      plot_image, bytscl(modelpsf[*,*,im], 0, 500)
      ;      stop
      skipimage:
    endfor  ; end loop over image
    
    ; remove zeros
    good=where(npix_psf gt 10, ngood, complement=bad)
    
    flag[bad]=0
    flag[good]=2
    
    n2=npix_psf[good]
    
    ;plot, n2, color=cgcolor('black'), psym=2
    res=n2-smooth(n2, 45, /edge_truncate)
    
    ;    plot, res, color=cgcolor('purple'), psym=2
    ;    stop
    rej=where(res le median(res)-(robust_sigma(res)*5) OR $
      res ge median(res)+(robust_sigma(res)*5), nrej, complement=keep)
      
;      plot, jd1[good], res, color=cgcolor('black'), psym=2
;      oplot, jd1[good[rej]], res[rej], color=cgcolor('red'), psym=2
;      stop
      
    ; images at xx are "good PSFs" - flag=2
    flag[good[keep]]=2
    flag[good[rej]]=0
    
    ; check flaged images
    ;    zeros=where(flag eq 0, nz)
    ;    for i=0, nz-1 do begin
    ;      plot_image, data1[*,*,zeros[i]], title=npix_psf[zeros[i]], color=cgcolor('black')
    ;      wait, 0.5
    ;    endfor
    ;    stop
    
    ;    twos=where(flag eq 2, nt)
    ;    for i=0, nt-1 do begin
    ;      plot_image, data1[*,*,twos[i]], title=npix_psf[twos[i]], color=cgcolor('black')
    ;      wait, 0.5
    ;    endfor
    ;stop
    
    ; save output
    fileout=outdir+file_basename(filein[f])
    
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag,  $
      psf_loc, npix_psf, modelpsf
      
  endfor  ; end loop over file
  
  
  print, 'end of program'
end