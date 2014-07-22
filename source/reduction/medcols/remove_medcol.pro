pro remove_medcol

; Do median column removal of p0 files using the full frame (medcol1) or raster (medcol2) values
;
Compile_opt idl2
  
sat='BA'
field='ORION'

opt='2' ; which median columns to use
  
indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'; CONTAINS ALL RASTER DIRECTORIES

tardir=file_search(indir+'HD*/', count=ntar)
  
outdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol'+opt+'/'

for tar=0, ntar-1 do begin  ; loop over each target directory
  filesin=file_search(tardir+'/*.sav', count=nf)  
  
  fname=file_basename(tardir[tar], '.sav') ; for saving output
  
  for i=0, nf-1 do begin  ; begin loop over each file
    
    ; restore the whole file - because we are re-saving
    restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
    
    if exp_time[0] ne 1. OR exp_ttl[0] ne 1. then continue ; do not process this file!
    
    nimg=n_elements(jd)
    
    tempdata=data1
    
    ; loop over each image and each column -
    ;   subtract the medcol1/medcol2 value, make any negative pixels non-zero
    for im=0, nimg-1 do begin
    
      data2=data1[*,*,im]
      
      xdim=(size(data2, /dim))[0]
      
      if opt eq '1' then for kk=0, xdim-1 do data2[kk,*]=data2[kk,*]-medcol1[kk,im]
      
      if opt eq '2' then for kk=0, xdim-1 do data2[kk,*]=data2[kk,*]-medcol2[kk,im]
      
      ; Check for negative pixels
      neg=where(data2 le 0, nneg)
      
      if nneg gt 0 then begin
      
        ineg=array_indices(data2, neg)
        
        ; replace negative pixels with a positive integer, not zero
        data2[ineg[0,*],ineg[1,*]]=1
        
      endif
      
      ; save new data2 in tempdata
      tempdata[*,*,im]=data2
      
    endfor  ; end loop over image - keep going until all images are processed
    
    ; rename tempdata to data1 and save result in new file
    data1=tempdata
    
    ; check that all exp_time an exp_ttl are 1. and 1.
    xx=where(exp_time ne 1. or exp_ttl ne 1., nxx)
    if nxx gt 0 then stop
    
    ; get dimensions of data from data
    sdata=size(data1, /dim)
    xdim=strtrim(sdata[0],2)
    ydim=strtrim(sdata[1],2)
    
    ; check output directory exists and make if not
    chk=file_search(outdir+fname, count=nchk)
    
    if nchk eq 0 then spawn, 'mkdir '+outdir+fname
    
    fileout=outdir+fname+'/'+fname+'_'+xdim+'_'+ydim+'_p1.sav'
    
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
      
  endfor  ; end loop over file
  
endfor  ; end loop over target directory



print, 'End of program'
print, 'Tackle HPs and CPs'

end