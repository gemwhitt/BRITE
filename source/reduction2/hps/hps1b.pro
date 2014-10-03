pro hps1b, sat, field, target

  ; Use the label region program to map obvious HPs above a threshold in the images - then clean image by image
  ;
  ; Input: p2 files
  ;
  ; Output: New p3 file AND HP map file (reduction/hp_map_labreg/)
  ;
  ; Modifications: modify flag according to medimg0 value
  ;    modify data1
  ;
  Compile_opt idl2
  
  ; sat='BA'
  ; field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p2/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/p3/'
    
  thr=100 ; threshold for finding bright regions
  
  ; check output directories exist and make accordingly
  chk1=file_search(outdir, count=n1)
  if n1 eq 0 then spawn, 'mkdir -p'+outdir
  
  filein=file_search(indir+target+'*.sav', count=nf)
    
  for ff=0, nf-1 do begin ; begin loop over files
  
    restore, filein[ff] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
    
    fname=file_basename(filein[ff],'.sav')
    
    ; update flag according to medimg0
    bad=where(medimg0 ge 5000, nbad)
    if nbad gt 0 then flag[bad]=0
    
    nimg=n_elements(jd)
    
    s=size(data1, /dim)
    
    ; define new arrays for cleaned data and hp maps
    hpmap=intarr(s[0],s[1],nimg)

    ; loop over images individually 
    for im=0, nimg-1 do begin
   ; print,  im
      ; first check flag ne 0
      if flag[im] eq 0 then continue
      
      ; extract image
      im0=data1[*,*,im]-50.
           
      ; check for negative pixels
      neg=where(im0 lt 0, nneg)
      if nneg gt 0 then begin
        neg2d=array_indices(im0, neg)
        im0[neg2d[0,*],neg2d[1,*]]=0
      endif
      
      ; add a border
      im1=intarr(s[0]+2,s[1]+2)
      im1[1:s[0],1:s[1]]=im0
      
      ; use label region to determine number of illuminated regions and number of PSF pixels
      r1=label_region(im1 ge thr, /all_neighbors)
      r1=r1[1:s[0],1:s[1]]                         ; trim off border
      
      ; Use histogram to check for number of regions (largest 1 entries are background and PSF)
      hist1=histogram(r1, binsize=1, locations=loc1, reverse_indices=ri1)
      
      if n_elements(hist1) eq 1 then begin  ; no target OR bad
        flag[im]=0
        goto, skipimage
      endif
      
      ; sort hist1 in decreasing order
      sort1=reverse(sort(hist1))
      
      ; check that there are enough pixels for a PSF
      if hist1[sort1[1]] lt 20 or hist1[sort1[1]] gt 120 then flag[im]=0
      if hist1[sort1[1]] lt 20 or hist1[sort1[1]] gt 120 then continue
      
      if n_elements(hist1) eq 2 then continue ; no HPs
      
      ; the number of results in hist1 are the number of regions (blobs) > thr, plus the background pixels
      ; so if n_elements(hist1) > 2 then - HP!
      if n_elements(hist1) gt 2 then begin
        
        map=intarr(s[0],s[1])
        
        ; record locations of HPs in hp map
        for region=2, n_elements(hist1)-1 do begin
          ; find 1d locations of HPs in each region
          ihp=ri1[ri1[sort1[region]]:ri1[sort1[region]+1]-1]
          
          nhp=n_elements(ihp)
          
          if nhp gt 3 then continue; not HPs - skip over
          
          ihp2d=array_indices(im0, ihp)
        
          xi=reform(ihp2d[0,*])
          yi=reform(ihp2d[1,*])
        
          ; update HP map
          map[xi,yi]=1
        
        endfor  ; end over region
        
     ;   plot_image, map
        
        hpmap[*,*,im]=map
;        stop
        
      endif ; if no HPs then nothing changed in image or added to HP map
      
      skipimage:
    endfor  ; end loop over images
    
    good=where(flag ne 0, ngood)
    hpmap=hpmap[*,*,good]
    
    ; split the hpmap into sections by time and check for frequency
    time1=jd[good]
    time2=time1[1:ngood-1]
    gap=where((time2-time1) gt 0.015, ngap)
    gap=[-1,gap,ngood-1]
    
    hp=[]
    hpmap2=data1*0
    
    for gp=0, ngap do begin
      
      iloc=indgen(gap[gp+1]-gap[gp])+(gap[gp]+1)
      ni=n_elements(iloc)
      
      if ni lt 10 then continue
      
      ; check frequency of HPs
      totmap=lonarr(s[0],s[1])
      
      for ii=0, ni-1 do totmap=totmap+hpmap[*,*,iloc[ii]]
      
      thr=ni*0.5
      
      hp=where(totmap ge thr, nhp)
      
      if nhp eq 0 then continue
      
      xi=(array_indices(totmap, hp))[0,*]
      yi=(array_indices(totmap, hp))[1,*]
      
      
      
; for im=0, ni-1 do begin
;        
;        plot_image, bytscl(data1[*,*,iloc[im]], 0, 500)
;        
;        oplot, xi, yi, color=cgcolor('purple'), psym=2
;      
;      wait, 2
;      
;   endfor
      
      map=totmap*0
      map[xi,yi]=1
      
      for jj=0, ni-1 do hpmap2[*,*,iloc[jj]]=map
      
endfor  ; end loop over group

; clean HPs
clean_hps, data1, hpmap2, flag

hpmap=hpmap2


    ; save p3 data
    fileout=outdir+fname+'_p3.sav'
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag, hpmap
      
      
  endfor  ; end loop over files
  
  print, 'End of program'
  print, 'Run hps2'
  
end