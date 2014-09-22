pro label_region_flag

; Program to determine frames with an adequate number of PSF pixels - use label_region.pro
; Use p2 images - which have been mostly - but not entirely cleaned of HPs
; 
; This modification of label_region_psf is for WERNER TESTSET files 
; 
; Output: Save a new variation of the input files, but in the date range of TOR and with flag added to denote
; good/bad/blank frames
;
Compile_opt idl2
  
indir='~/BRITE/TESTSETS/werner4lc/p2/'

filein=file_search(indir+'*.sav', count=nf)

outdir='~/BRITE/TESTSETS/werner4lc/p3/'
  
for f=0, 0 do begin ;nf-1 do begin
  
  print, file_basename(filein[f], '.sav')
    
  restore, filein[f]  ;jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag
  
  nimg=n_elements(jd)
  
  ; establish psf_loc array for recording the boundary of the PSF pixels - use to create a model?
  psf_loc=roi_loc*0
  npix_psf=intarr(nimg)
  ; hp_loc records any remaining HPs - in 2d as 0/1 (where 1=HP)
  hp_loc=data1*0.
  nhp=intarr(nimg)
  ; record exact PSF pixels and their values - save in separate array
  modelpsf=data1*0.
    
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
          
      ; the number of results in hist2 are the number of regions (blobs) > thr, plus the background pixels
      ; so if n_elements(hist2) > 2 then - HP!
      if n_elements(hist2) gt 2 then begin
        nhp[im]=total(hist2[sort2[2:n_elements(hist2)-1]])
        
        ; record locations of HPs
        for region=2, n_elements(hist2)-1 do begin
          ; find 1d locations of HPs in each section
          ihp=ri2[ri2[sort2[region]]:ri2[sort2[region]+1]-1]
          i2hp=array_indices(dat, ihp)
          for hp=0, n_elements(ihp)-1 do hp_loc[i2hp[0,hp], i2hp[1,hp], im]=1
;          plot_image, bytscl(dat, 0, 200)
;          oplot, i2hp[0,*], i2hp[1,*], color=cgcolor('orange'), psym=2
          
        endfor
      endif
     
      ; now overplot 2nd biggest group of pixels i.e. PSF pixels      
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
      for pix=0, npix_psf[im]-1 do modelpsf[xi[pix],yi[pix],im]=data1[xi[pix],yi[pix],im]
      
;      wset, 1
;      plot_image, bytscl(modelpsf[*,*,im], 0, 500)
;      stop
      skipimage:
    endfor  ; end loop over image
    
    ; remove zeros
    good=where(npix_psf gt 10, complement=bad)
    
    flag[bad]=0
    flag[good]=2
    
    n2=npix_psf[good]
    
    ;plot, n2, color=cgcolor('black'), psym=2
    res=n2-smooth(n2, 45, /edge_truncate)
    
;    plot, res, color=cgcolor('purple'), psym=2
;    stop
    xx=where(res ge median(res)-(robust_sigma(res)*3) AND $
      res le median(res)+(robust_sigma(res)*3), nxx, complement=yy)
    
    ; images at xx are "good PSFs" - flag=2  
    flag[good[xx]]=2
    flag[good[yy]]=0
    
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
    
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1, flag, $
    psf_loc, npix_psf, hp_loc, nhp, modelpsf
    
  endfor  ; end loop over file
  
  
  print, 'end of program'
end