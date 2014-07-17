pro save_medcol

; Calculate raster median columns, extract full-frame median columns
;   Re-save the data array - without the extra row
;   Save medcol1 (full-frame) and medcol2 (raster)
;   
;   Compile_opt idl2

sat='BA'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'

outdir=indir ; no point in keeping the original p0 files

filesin=file_search(indir+'*_p0*.sav', count=nfiles)  ; a or b - doesn't matter

for i=1, nfiles-1 do begin

  ; restore the while file
  restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
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
    
  endfor
  
  data1=tempdata
  
  ;plotsym, 0, /fill, 1.5
  ;!p.background=cgcolor('white')
  ;plot, ndead, color=cgcolor('black'), psym=8, xrange=[1500,1600], yrange=[0,1030], ystyle=1
  ;oplot, nsat, color=cgcolor('green'), psym=8
  ;stop
  
  ; save output
  fileout=filesin[i]
  ;stop
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
  
  
endfor

print, 'End of program'
print, 'Do remove_medcol1 OR remove_medcol2'

end
