pro hps2, sat, field, target

  ; Use the original find_hps method to find remainins HPs in the p3 images - close to or inside the PSF
  ; 
  ; This needs to be done on an image by image PLUS orbit by orbit basis
  ;
  ; Input: p3 files
  ;
  ; Output: New p4 file 
  ;
  ; Clean data1 by replacing HP and CPs with local median values
  ;
  Compile_opt idl2
  
 ; sat='BA'
 ; field='CENTAURUS'
  
  ;target=['HD127973','HD129056']
  ;target=['HD35411','HD37128']
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p3/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/p4/'
      
  ; check output directories exist and make accordingly
  chk1=file_search(outdir, count=n1)
  if n1 eq 0 then spawn, 'mkdir -p '+outdir
  
  filein=file_search(indir+target+'*.sav', count=nf)
  
  for ff=0, nf-1 do begin ; begin loop over files
  
    restore, filein[ff] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
       
    fname=file_basename(filein[ff],'_p3.sav')
   
    nimg=n_elements(jd)
        
    s=size(data1, /dim)
    
    ; define new arrays for cleaned data and hp maps
    newdata=data1
    
    ; add on a border to data1 
    data2=fltarr(s[0]+2,s[1]+2,nimg)
    data2[1:s[0],1:s[1],*]=data1
    
    ; work orbit by orbit
    jd1=jd-jd[0]
    jd2=jd1[1:nimg-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap)
    gap=[-1,gap,nimg-1]
    
    for orb=0, ngap do begin
      
      iloc=indgen(gap[orb+1]-gap[orb])+(gap[orb]+1)
      
      good=where(flag[iloc] eq 2, ngood, complement=bad)
      
      if ngood eq 0 then continue ; go to next orbit
      
      if n_elements(bad) gt ngood then begin
        flag[iloc]=0
        goto, skip_orbit
      endif
      
      iloc=iloc[good]
       
      hpmap=[]  ;define map array for this orbit
            
      ; loop over images in this orbit
      for im=0, ngood-1 do begin
              
        nhp=0
              
        im1=data1[*,*,iloc[im]]
        im2=data2[*,*,iloc[im]]
        
        ; check for negative pixels
        neg=where(im1 lt 0, nneg)
        if nneg gt 0 then begin
          neg2d=array_indices(im1, neg)
          im1[neg2d[0,*],neg2d[1,*]]=0
          
          neg=where(im2 lt 0, nneg)
          neg2d=array_indices(im2, neg)
          im2[neg2d[0,*],neg2d[1,*]]=0
        endif
        
        ; wset, 0
        ; plot_image, bytscl(im0, 0, 500)
        ;
        ; identify ALL "bright" pixels, i.e. DN > thr 
        thr=50
        
        br=where(im2 gt thr, nbr)
        
        if nbr eq 0 then begin
          flag[iloc[im]]=0
          goto, skipimage
        endif
        
        ; find 2d coords of bright pixels in the padded image
        xi=(array_indices(im2, br))[0,*]
        yi=(array_indices(im2, br))[1,*]
        
        ; criteria 2
        cr2=3
        
        map1=fltarr(s[0],s[1])  ; original dimensions
        
        for p=0, nbr-1 do begin
        
          x1=(xi[p]-1) 
          x2=(xi[p]+1) 
          y1=(yi[p]-1) 
          y2=(yi[p]+1) 
          
          cutout = im2[x1:x2,y1:y2]     ; 3x3 neighbourhood
          
          ; calculate median of cutout
          localmed=median(cutout)
          
          j = where((cutout-localmed) ge 3*robust_sigma(cutout)+localmed, nj)       ; by cr2 above local median
          
          if nj eq 0 then continue  ; no HPs here
          
          cj=where(j eq 4)  ; check that central pixel in square is the highest
          
          ; this is NOT a HP if...
          if cj[0] eq -1 OR nj gt cr2 then continue
          
;          plot_image, bytscl(im2, 0, 500), title=iloc[im], charsize=0.9, color=cgcolor('black')
;          oplot, [xi[p]], [yi[p]], color=cgcolor('purple'), psym=2
;          stop
          
          ; modify xi and yi
          xi[p]=xi[p]-1
          yi[p]=yi[p]-1
          
          ; otherwise, record HP location in hpmap  
          map1[xi[p],yi[p]]=im1[xi[p],yi[p]]
          
          nhp=nhp+1           
          
        endfor  ; end loop over bright pixels
        
        ; add map1 to hpmap
        hpmap=[[[hpmap]],[[map1]]]
        
;        wset, 0
;        plot_image, bytscl(hpmap[*,*,iloc[im]], 0, 200)
;        
;        stop
             
        skipimage:
      endfor  ; end loop over images
      
      ; check for consistency of HPs... must be in > 50% of images to be a HP
      
      ; add images together and check for number of potential HPs
      if n_elements(hpmap) eq 0 then continue
      
      temp_map=im1*0
      smap=size(hpmap, /dim)
      if n_elements(smap) eq 2 then temp_map=hpmap else for ii=0, smap[2]-1 do temp_map=temp_map+hpmap[*,*,ii]
      
      xx=where(temp_map gt 0, nhp)
      
  ;    plot_image, bytscl(temp_map, 0, 500)
        
      
      if nhp eq 0 then continue
      
      for hp=0, nhp-1 do begin
        
        ; get 2d location of this pixel in hpmap
        xi=(array_indices(temp_map, xx[hp]))[0]
        yi=(array_indices(temp_map, xx[hp]))[1]
        
        freq=where(hpmap[xi,yi,*] gt 0, nfreq)
        
        if n_elements(smap) eq 2 then pcent=float(nfreq)*100. else pcent=float(nfreq)/float(smap[2])*100.
            
        if pcent ge 20. then begin
        
;        plot_image, bytscl(temp_map, 0, 500), title=iloc[0], color=cgcolor('black'), charsize=0.7
;        oplot, [xi], [yi], color=cgcolor('purple'), psym=2
;        print, pcent
;        stop

          ; clean images in this orbit of HP and CP in each image 
          for ii=0, ngood-1 do begin
            
            im1=data1[*,*,iloc[ii]]
            
            x1=xi-1 > 0
            x2=xi+1 < s[0]-1
            y1=yi-1 > 0
            y2=yi+1 < s[0]-1
            
            cutout=im1[x1:x2,y1:y2]
            
            im1[xi,yi]=median(cutout) ; this cleans the HP
            
            if xi lt s[0]-2 then im1[xi+1,yi]=median(im1[xi:xi+2,y1:y2])
            
            data1[*,*,iloc[ii]]=im1 

          endfor
                   
        endif

        
      endfor
      skip_orbit:
    endfor  ; end loop over orbit
    
;    for im=0, nimg-1, 100 do begin
;      if flag[im] ne 2 then continue
;      
;      plot_image, bytscl(data1[*,*,im], 0, 500)
;      wait, 2
;    endfor

    ; save p4 data
    fileout=outdir+fname+'_p4.sav'
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
      
      
  endfor  ; end loop over files
  
  print, 'End of program'
  print, 'do get_psf_boundary and get_com'
  
end