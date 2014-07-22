pro save_medcol

;   Calculate raster median columns, extract full-frame median columns
;   Re-save the data array - without the extra row
;   Save medcol1 (full-frame) and medcol2 (raster)
;   
Compile_opt idl2

sat='BA'
field='CENTAURUS'

indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'

tardir=file_search(indir+'HD*/', count=ntar) ; target directories

for tar=0, ntar-1 do begin
  
  filesin=file_search(tardir[tar]+'/*.sav', count=nf) 
  
  for ff=0, nf-1 do begin
  
    ; restore the whole file
    restore, filesin[ff]   ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype
    
    nimg=n_elements(jd)
    
    xdim=(size(data1, /dim))[0]
    ydim=(size(data1, /dim))[1]
    
    medcol1=lonarr(xdim, nimg)  ; array of full-frame median columns - calculated onboard
    medcol2=lonarr(xdim, nimg)  ; array of raster-calculated column medians
    
    medimg0=fltarr(nimg)        ; array of whole image medians (raw)
    
    ndead=intarr(nimg)  ; initial number of dead pixels on image
    nsat=intarr(nimg)   ; number of pixels above a given threshold on the image
    
    tempdata=lonarr(xdim, ydim-1,nimg) ; temporary array for new dataset with extra row trimmed off
    
    for im=0, nimg-1 do begin
    
      data2=data1[0:xdim-1,0:ydim-2,im]   ; image data
      
      tempdata[*,*,im]=data2
      
      medimg0[im]=median(data2)           ; calculate image median
      
      dead=where(data2 le 0, num_dead)    ; check for dead pixels
      ndead[im]=num_dead
      
      sat=where(data2 ge 10000, num_sat)  ; check for nonlinear
      nsat[im]=num_sat
      
      ; extract the fullframe median columns
      medcol1[*,im]=data1[*,ydim-1,im]
      
      ; calculate median from raster
      for k=0, xdim-1 do medcol2[k,im]=median(data2[k,0:ydim-2])
      
    endfor  ; end loop over this image
      
    data1=tempdata
      
    ; save output
    fileout=filesin[ff]
  
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
   
  
endfor  ; end loop over this file
  
endfor  ; end loop over this target

print, 'End of program'
print, 'Do remove_medcol1 OR remove_medcol2'

end
