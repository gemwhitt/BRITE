pro map_hps_test

  ; PURPOSE - same as map_hps1 - but for testing
  ; Save both hps1 and hps2 along with data1 - then analyse the results
  ; Output: hps1 and hps2 make Hot pixel and CP maps for each image
  
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 1.3
  
  sat='BA'
  field='ORION'
  thr=50  ; threshold above which HPs are identified

  indir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/'
  filesin=file_search(indir+'*.sav', count=nf)
  
  outdir='~/BRITE/'+sat+'/'+field+'/reduction/hp_maps/'
  
  for f=0, 0 do begin ;nf-1 do begin
  
    obj=obj_new('IDL_Savefile', filesin[f])
    obj->restore, 'jd'
    obj->restore, 'data1'
    obj->restore, 'ccd_temp'
    obj->restore, 'vmag'
    obj->restore, 'medimg0'
    
    jd1=jd-jd[0]
    
    nfrm=n_elements(jd)
    
    fname=file_basename(filesin[f],'.sav')
    
    totdn=lonarr(nfrm)
    for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
    
    sdata=size(data1, /dim)
    xdim=sdata[0]
    ydim=sdata[1]
    
    hps1=lonarr(xdim, ydim, nfrm)
    hps2=lonarr(xdim, ydim, nfrm)
    flag=intarr(nfrm)+2 ; 0=bad, 1=blank, 2=target - to begin assume everything is target
    
    jd2=jd1[1:nfrm-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits
    gap=[-1,gap,nfrm-1]
    
    
    for orbit=0, 50 do begin  ;ngap do begin                                  ; LOOP OVER EACH ORBIT
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1         ; INDICES OF IMAGES IN THIS ORBIT
      
      imgs=data1[*,*,iloc]
      
      nimg=n_elements(iloc)
      
      rej1=where(medimg0[iloc] gt 5000 OR totdn[iloc] lt 5000, nrej1, complement=keep1) ; bad images - discard
      
      ngood=n_elements(keep1)
      
      if nrej1 gt 0 then begin                                  ;stop - analyse discarded images
      
        ; change flag for this image to 0=bad
        flag[iloc[rej1]]=0
        
        if ngood lt 2 then begin    ; IF there are less than 2 good frames then save previous HP values and go to next orbit
          for p=0, ngood-1 do hps[*,*,iloc[keep1[p]]]=hps[*,*,iloc[0]-1]
          
          goto, next_orbit
        endif
        
        ; OTHERWISE remove bad images from the analysis and continue
        iloc=iloc[keep1]
      endif
      
      ; Start Original method
      hpo=lonarr(xdim, ydim, ngood)                         ; array to collect HPs from each frame in this orbit     
      for img=0, ngood-1 do begin
      
        im=imgs[*,*,keep1[img]]
        
        s=size(im, /dim) ; length of array in x and y
          cr1=thr
          cr2=3
          ima=find_hps2(im,cr1,cr2,iloc[keep1[img]], ww,wx,wy)
          ;plot_image, bytscl(im, 0, 200)
          ;oplot, [wx], [wy], color=cgcolor('orange'), psym=2
          
          nhp=n_elements(wx)
          
          if nhp eq 0 then stop
          
          if nhp gt 0 then begin
            for hp=0, nhp-1 do begin
              hpo[wx[hp],wy[hp],img]=im[wx[hp],wy[hp]]
              ; check that CP is not off the image and record if available
              if wx[hp] le s[0]-2 then hpo[wx[hp]+1,wy[hp],img]=im[wx[hp]+1,wy[hp]]
            endfor
          endif
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;       
      endfor  ; end loop over this image
      
      ; the result from either option is hpo - array with HPs and CPs saved for all images in this orbit
      
      ; check frequency of HP occurance for this orbit
      for ii=0, xdim-1 do begin  ; check every pixel in x and y
        for jj=0, ydim-1 do begin
          pix=hpo[ii,jj,*]
          
          xx=where(pix gt thr, nxx)
          
          if nxx ge 0.25*ngood then begin $ ; HP!
            hps1[ii,jj,iloc]=hpo[ii,jj,*]
          if ii le xdim-2 then hps1[ii+1,jj,iloc]=hpo[ii+1,jj,*]
        endif
      endfor
    endfor
    
    ; Start new method
    hpo=lonarr(xdim, ydim, ngood)                         ; array to collect HPs from each frame in this orbit
    
    for img=0, ngood-1 do begin
    
      im=imgs[*,*,keep1[img]]
      
      s=size(im, /dim) ; length of array in x and y

    ; NEW METHOD....
    im1=lonarr(s+2)   ; append array with extra row/column (for label_region.pro)
    im1[1,1]=im      ; ""
    
    r1=label_region(im1 gt thr, /all_neighbors)  ; find groups of pixels above thr
    r1=r1[1:s[0],1:s[1]]                         ; trim off border
    
    ; do a histogram of r1 and save indices in ri
    result=histogram(r1, binsize=1, reverse_indices=ri, locations=loc)
    
    ; ignore any groups of pixels with > than 3 pixels
    keep=where(result le 3, nkeep, complement=rej)    ; this should ignore the background pixels and the target pixels
    
    IF nkeep eq 0 then stop ;CONTINUE ; go to next image
    
    i=0
    repeat begin
      ind=ri[ri[keep[i]]:ri[keep[i+1]-1]] ; 1d index of HPs - convert to 2d
      
      ind2=array_indices(r1, ind)
      
      hpo[ind2[0],ind2[1],img]=im[ind2[0],ind2[1]]
      
      ; check that CP is not off the image and record if available
      if ind2[0] le s[0]-2 then hpo[ind2[0]+1,ind2[1],img]=im[ind2[0]+1,ind2[1]]
      
      i=i+1
    endrep until i eq nkeep-1
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
    endfor
    
    ; check frequency of HP occurance for this orbit
    for ii=0, xdim-1 do begin  ; check every pixel in x and y
      for jj=0, ydim-1 do begin
        pix=hpo[ii,jj,*]
        
        xx=where(pix gt thr, nxx)
        
        if nxx ge 0.25*ngood then begin $ ; HP!
          hps2[ii,jj,iloc]=hpo[ii,jj,*]
        if ii le xdim-2 then hps2[ii+1,jj,iloc]=hpo[ii+1,jj,*]
      endif
    endfor
  endfor
    
  next_orbit:
    
  endfor ; end loop over this orbit
 
  fileout=outdir+fname+'_hpmaps.sav'
  
  save, filename=fileout, data1, hps1, hps2
  
endfor  ;end loop over this file

print, 'end of program'
end