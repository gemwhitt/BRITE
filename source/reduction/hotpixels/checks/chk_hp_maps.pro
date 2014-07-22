pro chk_hp_maps

; Remove HPs from the images and cycle though to check for "goodness" of identification.
; 
; Run after map_hps_test.pro
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.3

sat='BA'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/reduction/hp_maps/'
filesin=file_search(indir+'*hpmaps.sav', count=nf)

for f=0, 0 do begin ;nf-1 do begin
  
  restore, filesin[f] ;data1, hps1, hps2
  
  ; redefine data1 as data0 - to free up data1
  data0=data1
  
  ; loop over each image in data0 and make cleaned images in data1 and data2 using hps1 and hps2
  data1=data0*0
  data2=data0*0
  
  sdata=size(data0, /dim)
  
  nimg=sdata[2]
  
  for im=0, nimg-1 do begin
    
      
    ; DO ORIGINAL METHOD FIRST  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    im1=data0[*,*,im]
    hp1=hps1[*,*,im]
    
    ind=where(hp1 gt 0, nxy)
    if nxy eq 0 then continue
    
    xy=array_indices(hp1, ind)
    
    if nxy eq 0 then continue
    
    xi=xy[0,*]
    yi=xy[1,*]
    
    for ii=0, nxy-1 do begin
      x1=xi[ii]-1 > 0
      x2=xi[ii]+1 < (sdata[0]-1)
      y1=yi[ii]-1 > 0
      y2=yi[ii]+1 < (sdata[1]-1)
      
      im1[xi[ii],yi[ii]]=median(im1[x1:x2,y1:y2])
      
      ; and do CPs
      if xi[ii] le (sdata[0]-3) then im1[xi[ii]+1,yi[ii]]=median(im1[x1+1:x2+1,y1:y2])
    endfor
    
    data1[*,*,im]=im1
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    ; NOW DO NEW METHOD ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    im2=data0[*,*,im]
    hp2=hps2[*,*,im]
    
    ind=where(hp2 gt 0, nxy)
    if nxy eq 0 then continue
    
    xy=array_indices(hp2, ind)
    
    xi=xy[0,*]
    yi=xy[1,*]
    
    for ii=0, nxy-1 do begin
      x1=xi[ii]-1 > 0
      x2=xi[ii]+1 < (sdata[0]-1)
      y1=yi[ii]-1 > 0
      y2=yi[ii]+1 < (sdata[1]-1)
      
      im2[xi[ii],yi[ii]]=median(im2[x1:x2,y1:y2])
      
      ; and do CPs
      if xi[ii] le (sdata[0]-3) then im2[xi[ii]+1,yi[ii]]=median(im2[x1+1:x2+1,y1:y2])
    endfor
    
    data2[*,*,im]=im2  
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
  endfor  ; end loop over images
  
  ; now plot both
  for im=0, nimg-1 do begin
    
    
    wset, 0
    plot_image, bytscl(data0[*,*,im], 20, 500)
    
    wset, 1
    plot_image, bytscl(data1[*,*,im], 20, 500)
    
    wset, 2
    plot_image, bytscl(data2[*,*,im], 20, 500)
    
stop
  endfor
  
  
  stop
  
  
  
  
  
endfor




stop
end

