pro remove_medcol

; Do median column removal of p0 files using the full frame (medcol1) or raster (medcol2) values
;
Compile_opt idl2
  
sat='BA'
field='ORION'

opt='2'
  
indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'
  
outdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol'+opt+'/'
  
filesin=file_search(indir+'*_p0.sav', count=nfiles)  ; a or b - doesn't matter
  
for i=0, nfiles-1 do begin
  
  fname=file_basename(filesin[i], '_p0.sav') ; change this for p0b
  
  ; restore the while file
  restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
                        ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
    
  nimg=n_elements(jd)
  
  tempdata=data1
  
  ; loop over each image and each column - 
  ;   subtract the medcol1 value, i.e. the full frame value, make any negative pixels non-zero
  for im=0, nimg-1 do begin
    
    data2=data1[*,*,im]
    
    xdim=(size(data2, /dim))[0]
    
    if opt eq '1' then for kk=0, xdim-1 do data2[kk,*]=data2[kk,*]-medcol1[kk,im]
    
    if opt eq '2' then for kk=0, xdim-1 do data2[kk,*]=data2[kk,*]-medcol2[kk,im]
    
    ; Check for negative pixels
    neg=where(data2 le 0, nneg)
    
    if nneg gt 0 then begin
      
      ineg=array_indices(data2, neg)
      
      ; replace negative pixels with a positive integer
      data2[ineg[0,*],ineg[1,*]]=1
      
    endif
    
    ; save new data2 in tempdata
    tempdata[*,*,im]=data2    
    
  endfor
  
  ; rename data2 to data1 and save result in new file - appended with _p1a or _p1b
  data1=data2
  
  fileout=outdir+fname+'_p1.sav'
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
  
endfor

print, 'End of program'
print, 'Tackle HPs and CPs'

end