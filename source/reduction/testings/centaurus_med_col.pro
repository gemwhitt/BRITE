pro centaurus_med_col

  ; Program to remove warm columns from the rasters - values for column medians are saved in the data - at the end of columns
  ;
  Compile_opt idl2
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  sat='UB'
  
  indir='~/BRITE/data/'+sat+'/roi_raw_sav/CENTAURUS/'
  
  outdir='~/BRITE/data/UB/p1/CENTAURUS/'; location where save files will go
  
  tar=['HD122980', 'HD138690']
  
  filesin=file_search(indir+'*'+tar+'*.sav', count=nfiles)
  
  for i=1, nfiles-1 do begin
  
    fname=file_basename(filesin[i], '_p0.sav')
    
    restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl
    
    data2=data1[0:27,0:28,*]
    data1=data2
    
    nfrm=n_elements(jd)
    
    medcol1=fltarr(28,nfrm)      ; medcol values from the whole image
    medcol2=fltarr(28,nfrm)      ; median value of the raster column
    
    for j=0, 500 do medcol1[*,j]=data2[*,28,j]
    for j=0, 500 do for k=0, 27 do medcol2[k,j]=median(reform(data2[k,0:27,j]))
    
    plotsym, 0, /fill, 0.8
 ;   window, 0, xsize=600, ysize=550
    
  ;  for j=0, 63 do begin
      
  ;    print, j

     ; if j eq 0 then begin
     ;   plot, medcol1[*,j], medcol2[*,j], color=cgcolor('black'), psym=8, xrange=[0,500], yrange=[0,500] 
     ;   oplot, [0,500], [0,500], thick=3, color=cgcolor('blue')
     ; endif else begin
     ;   if j lt 40 then oplot, medcol1[*,j], medcol2[*,j], color=cgcolor('black'), psym=8 else $
     ;     oplot, medcol1[*,j], medcol2[*,j], color=cgcolor('red'), psym=8
     ;  endelse
 
   ; endfor


  ;  stop
    for j=50, 50 do begin
      wset, 0
      plot_image, bytscl(data1[*,*,j], 20, 1000)
      wait, 2
    endfor
    
    print, median(data1[*,0:27,50])
    print, robust_sigma(data1[*,0:27,50])
    
    data3=reform(data1[*,0:27,50])
    for k=0, 27 do data3[k,*]=data3[k,*]-medcol1[*,50]
    
    wset, 1;window, 1, xsize=600, ysize=550
    plot_image, bytscl(data3, 20, 1000)
    
    print, median(data3)
    print, robust_sigma(data3)
    
    
    
    data3=reform(data1[*,0:27,50])
    for k=0, 27 do data3[k,*]=data3[k,*]-medcol2[*,50]
    
    wset, 2 ;window, 2, xsize=600, ysize=550
    plot_image, bytscl(data3, 20, 1000)
    
    print, median(data3)
    print, robust_sigma(data3)
    
    
    stop
    for j=50, 50 do begin
      plot_image, bytscl(data1[*,*,j], 20, 1000)
      wait, 2
    endfor
    stop
    
    nimg=(size(data1, /dim))[2]
    
    xdim=(size(data1, /dim))[0]
    ydim=(size(data1, /dim))[1]
    
    medcols=fltarr(xdim-1, nimg)
    
    medimg=fltarr(nimg)
    
    rightcol=fltarr(xdim-1,nimg)
    
    tempdata=lonarr(xdim-1, ydim-1,nimg)
    
    for j=0, nimg-1 do begin
    
      data2=data1[0:xdim-2,0:ydim-2,j]
      stop
      ; check upper left corner for a non-zero value
      chkpix=data1[0,ydim-1,j]
      
      if chkpix gt 0 then medcols[*,j]=data1[0:xdim-2,ydim-1,j] else begin
      
        for k=0, xdim-2 do medcols[k,j]=median(data1[k,0:ydim-2])
        
      endelse
      
      ; calculate image median
      medimg[j]=median(data2)
      
      ; check for dead pixels
      dead=where(data2 le 0.0, ndead)
      
      ; check if there is data in the righthand column
      rightcol[*,j]=reform(data1[xdim-1,0:ydim-2,j])
      
      for k=0, xdim-2 do tempdata[k,*,j]=data2[k,*]-medcols[k,j]
      
      ; make any negative values zero values instead
      data2=reform(tempdata[*,*,j])
      
      xx=where(data2 lt 0, nneg)
      
      if nneg gt 0 then begin
      
        loc2d=array_indices(data2, xx)
        
        data2[loc2d[0,*], loc2d[1,*]]=0.0
        
        xx=where(data2 lt 0, nneg)
        
        if nneg gt 0 then stop
        
        tempdata[*,*,j]=data2
        
      endif
      
    endfor
    
    data1=tempdata
    
    ; save data
    fileout=outdir+fname+'_p1.sav'
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
      medcols, medimg, ndead, rightcol
      
      
      
  endfor
  
  
  print, 'end of program'
end


