pro hps1

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
  
  sat='UB'
  field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p2/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/p3/'
  
  mapdir='~/BRITE/'+sat+'/'+field+'/reduction/hp_maps/'
  
  thr=50  ; threshold for finding bright regions
  
  ; check output directories exist and make accordingly
  chk1=file_search(outdir, count=n1)
  if n1 eq 0 then spawn, 'mkdir -p'+outdir
  
  chk2=file_search(outdir, count=n2)
  if n2 eq 0 then spawn, 'mkdir -p'+mapdir
  
  filein=file_search(indir+'*.sav', count=nf)
  
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
    hp_map=data1*0
    newdata1=data1
    nhp1=intarr(nimg)
    
    ; add border onto each image for label_region program to work
    data2=fltarr(s[0]+2,s[1]+2,nimg)
    data2[1:s[0],1:s[1],*]=data1
    
    ; loop over images individually (unless flag[im]=0), record HPs in hp_map arrays - then clean image - save new image
    for im=0, nimg-1 do begin
    
      ; first check flag ne 0
      if flag[im] eq 0 then continue
      
      ; extract image
      im0=data1[*,*,im]
      im1=data2[*,*,im]
      map=im0*0
 
      ; check for negative pixels
      neg=where(im0 lt 0, nneg)
      if nneg gt 0 then stop
      
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
      
      ; record median of background pixels
      bkgd_ind=ri1[ri1[sort1[0]]:ri1[sort1[0]+1]-1]
      bkgd_2di=array_indices(im0, bkgd_ind)
      bkgd=median(data1[bkgd_2di[0,*],bkgd_2di[0,*],im]) > 0
      
      ; the number of results in hist1 are the number of regions (blobs) > thr, plus the background pixels
      ; so if n_elements(hist1) > 2 then - HP!
      if n_elements(hist1) gt 2 then begin
        nhp1[im]=total(hist1[sort1[2:n_elements(hist1)-1]])
        
        ihp=[]
        
        ; record locations of HPs in hp map
        for region=2, n_elements(hist1)-1 do $
          ; find 1d locations of HPs in each region
          ihp=[ihp,ri1[ri1[sort1[region]]:ri1[sort1[region]+1]-1]]
          ihp_2d=array_indices(im0, ihp)
          
          ; check hp_map versus image
;          wset, 0
;          plot_image, bytscl(im0, 0, 500)
;          plotsym, 0, /fill, 1.2
;          oplot, ihp_2d[0,*], ihp_2d[1,*], color=cgcolor('purple'), psym=8
;          
;          wset, 1
;          plot_image, bytscl(hp_map[*,*,im], 0, 500)
          
            xi=reform(ihp_2d[0,*])
            yi=reform(ihp_2d[1,*])
            
            ; update HP map
            map[xi, yi]=im0[xi, yi]
            hp_map[*,*,im]=map
            
;            wset, 1
;            plot_image, bytscl(hp_map[*,*,im], 0, 500)
;            stop
            ; replace HP with bkgd
            im0[xi,yi]=bkgd
            newdata1[*,*,im]=im0                              
        
         ;check before and after image
;        wset, 0
;        plot_image, bytscl(data1[*,*,im], 0, 500)
;        
;        wset, 1
;        plot_image, bytscl(newdata1[*,*,im], 0, 500) 
;        stop
    
      endif ; if no HPs then nothing changed in image or added to HP map
     
      skipimage:
    endfor  ; end loop over images
    
    ; resave newdata1 as data1
    data1=newdata1
    
    ; save hp_map 
    mapout=mapdir+fname+'_hpmap1.sav'
    save, filename=mapout, hp_map
    
    ; save p3 data
    fileout=outdir+fname+'_p3.sav'
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag, nhp1
    
    
  endfor  ; end loop over files
  
  print, 'End of program'
  print, 'Run hps2'
  
end