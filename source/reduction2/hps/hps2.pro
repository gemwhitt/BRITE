pro hps2

  ; Use the original find_hps method to find remainins HPs in the p3 images - close to or inside the PSF
  ;
  ; Input: p3 files
  ;
  ; Output: New p4 file 
  ;
  ; Clean data1 by replacing HP and CPs with local median values
  ;
  Compile_opt idl2
  
  sat='UB'
  field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p3/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/p4/'
      
  ; check output directories exist and make accordingly
  chk1=file_search(outdir, count=n1)
  if n1 eq 0 then spawn, 'mkdir -p '+outdir
  
  filein=file_search(indir+'*.sav', count=nf)
  
  for ff=0, nf-1 do begin ; begin loop over files
  
    restore, filein[ff] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
    
    fname=file_basename(filein[ff],'.sav')
   
    nimg=n_elements(jd)
    
    s=size(data1, /dim)
    
    ; define new arrays for cleaned data and hp maps
    newdata1=data1
    
    ; add on a border to data1 
    data2=fltarr(s[0]+2,s[1]+2,nimg)
    data2[1:s[0],1:s[1],*]=data1
    
    ; loop over images, ignore images with flag=0, if there are no good images then skip over
    for im=0, nimg-1 do begin
      
      if flag[im] eq 0 then continue
      
      im0=data2[*,*,im]

      ; check for negative pixels
      neg=where(im0 lt 0, nneg)
      if nneg gt 0 then stop
         
;      wset, 0
;      plot_image, bytscl(im0, 0, 500)
;        
      ; identify ALL pixel above thr in the averaged image
      thr=50
      br=where(im0 gt thr, nbr)
      
      if nbr eq 0 then begin
        flag[im]=0
        goto, skipimage
      endif
      
      xi=reform((array_indices(im0, br))[0,*])
      yi=reform((array_indices(im0, br))[1,*])
      
      cr2=2
            
      for p=0, nbr-1 do begin
        
        x1=(xi[p]-1) 
        x2=(xi[p]+1) 
        y1=(yi[p]-1)
        y2=(yi[p]+1) 
        
        
        cutout = im0[x1:x2,y1:y2]     ; 3x3 neighbourhood
        
        ; calculate median of cutout
        localmed=median(cutout)
        
        j = where((cutout-localmed) gt 5*robust_sigma(cutout)+localmed, nj)       ; by cr2 above local median
        
        if nj eq 0 then continue  ; no HPs here
        
        cj=where(j eq 4)  ; check that central pixel in square is the highest
        
        ; this is NOT a HP if...
        if cj[0] eq -1 OR nj gt cr2 then continue
        
        ; otherwise, record HP location
;        plotsym, 0, /fill, 1.2
;        oplot, [xi[p]], [yi[p]], color=cgcolor('green'), psym=8
;        
        ; replace hp at xi, yi with median of cutout, then replace at x1+1 (if possible
        im0[xi[p],yi[p]]=localmed
        
        if xi[p] lt s[0]-3 then im0[xi[p]+1,yi[p]]=median(im0[xi[p]:xi[p]+2,yi[p]-1:yi[p]+1])

      endfor  ; end loop over bright pixels
      
;      wset, 1
;      plot_image, bytscl(im0, 0, 500)
;      stop

      ; trim border off im0
      im0=im0[1:s[0],1:s[1]]
      
      newdata1[*,*,im]=im0
          
        
        ; check before and after image
;                wset, 0
;                plot_image, bytscl(data1[*,*,im], 0, 500)
;        
;                wset, 1
;                plot_image, bytscl(newdata1[*,*,im], 0, 500)
;                stop
        
      
      skipimage:
    endfor  ; end loop over orbit
    
    ; resave newdata1 as data1
    data1=newdata1
    
    ; save p4 data
    fileout=outdir+fname+'_p4.sav'
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag, nhp1
      
      
  endfor  ; end loop over files
  
  print, 'End of program'
  print, 'do label region to get PSF boundary and COM'
  
end