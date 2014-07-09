pro median_column2

; PURPOSE: calculate raster median columns - save as medcol1 array [xdim,nfrm]
;          also extract full-frame median columns - save as medcol2 array [xdim,nfrm]
;          calculate median image p0 - save as medimg0[nfrm]
;          save med-col removed data as data1 & data2 - depending on method of removal
;          
Compile_opt idl2
  
sat='BA'
field='ORION'
  
indir='~/BRITE/data/'+sat+'/'+field+'/roi_raw/sav/'
  
outdir='~/BRITE/data/'+sat+'/'+field+'/p1/sav/'
  
filesin=file_search(indir+'*_p0b.sav', count=nfiles)
  
for i=0, nfiles-1 do begin
  
  fname=file_basename(filesin[i], '_p0b.sav')
  print, fname
    
  restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
                        ;simbad_radec, vmag, bmag, parlax, otype, sptype
  
  data0=data1
    
  nimg=n_elements(jd)
    
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
    
  medcol1=lonarr(xdim, nimg)
 ; medcol2=lonarr(xdim, nimg)
    
  medimg0=fltarr(nimg)

  ndead=intarr(nimg)  ; initial number of dead pixels on image
  nsat=intarr(nimg)   ; number of pixels above a given threshold on the image
    
  tempdata1=lonarr(xdim, ydim-1,nimg)
 ; tempdata2=lonarr(xdim, ydim-1,nimg)
  
  ; FOR DISPLAY  
  ; window,0, xsize=600, ysize=550, xpos=1500, ypos=100
  ; for im=0, nfrm-1 do begin
  ; plot_image, bytscl(data1[*,*,im], 20, 500)
  ; stop / wait, 0.5
  ; endfor
  ;stop
  
  for im=0, nimg-1 do begin
    
    data2=data1[0:xdim-1,0:ydim-2,im]   ; image data 
        
    medimg0[im]=median(data2)           ; calculate image median
     
    dead=where(data2 le 0, num_dead)    ; check for dead pixels
    ndead[im]=num_dead
    
    sat=where(data2 ge 10000, num_sat)  ; check for nonlinear
    nsat[im]=num_sat
      
    ; calculate median from raster 
    for k=0, xdim-1 do medcol1[k,im]=median(data2[k,0:ydim-2])
    
    ; save median value from full frame
  ;  medcol2[*,im]=data1[*,ydim-1,im]
      
    ; subtract the median columns
    for k=0, xdim-1 do tempdata1[k,*,im]=data2[k,*]-medcol1[k,im] ; raster median removed
    
;    for k=0, xdim-1 do tempdata2[k,*,im]=data2[k,*]-medcol2[k,im] ; full frm median removed
  
  endfor 
  
  data1=tempdata1
;  data2=tempdata2
  
  ; save output
  fileout=outdir+fname+'_p1b.sav'
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, data1
  
endfor
  
print, 'End of program'
print, 'Do HP mitigation'

end


