pro remove_medcol_1s

; Do median column removal of p0 files using the full frame (medcol1) or raster (medcol2) values
;
; Modified from remove_medcol.pro - process ONLY 1s files
; Save as: e.g. HDxxxx_24_24_p1_0.sav 
; 
; Remove top row from rasters which contain the full-frame  med cols
Compile_opt idl2
  
sat='BA'
field='CENTAURUS'
  
opt='2' ; which median columns to use
 
indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/2014_0607/'; CONTAINS ALL RASTER DIRECTORIES
  
tardir=file_search(indir+'HD*/', count=ntar)
  
outdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol'+opt+'/2014_0607/'
  
for tar=0, ntar-1 do begin  ; loop over each target directory
  
    tar_name=file_basename(tardir[tar], '.sav') ; name of target
  
    filesin=file_search(tardir[tar]+'/*.sav', count=nf)
    
    fname=file_basename(filesin, '.sav')
    
    for i=0, nf-1 do begin  ; begin loop over each file
    
      ; restore the whole file - because we are re-saving
      restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
      ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
      
      ; get index of this file from fname
      flen=strlen(fname[i])
      ind=strmid(fname[i], flen-1)
      
      if exp_time[0] ne 1. OR exp_ttl[0] ne 1. then continue ; do not process this file!
      
      nimg=n_elements(jd)
      
      tempdata=[]
      
      ; loop over each image and each column -
      ;   subtract the medcol1/medcol2 value, make any negative pixels non-zero
      for im=0, nimg-1 do begin
      
        data2=data1[*,*,im]
        
        xdim=(size(data2, /dim))[0]
        ydim=(size(data2, /dim))[1]
        
        if ydim mod 2 eq 1 then data2=data2[*,0:ydim-2] ; this removes the top row (full frame medcols)
        ;stop ; check this
        
        if opt eq '1' then for kk=0, xdim-1 do data2[kk,*]=data2[kk,*]-medcol1[kk,im]
        
        if opt eq '2' then for kk=0, xdim-1 do data2[kk,*]=data2[kk,*]-medcol2[kk,im]
        
        ; Check for negative pixels
        neg=where(data2 lt 0, nneg)
        
        if nneg gt 0 then begin
        
          ineg=array_indices(data2, neg)
          
          ; replace negative pixels with a positive integer, not zero
          data2[ineg[0,*],ineg[1,*]]=1
          
        endif
        
        ; save new data2 in tempdata
        tempdata=[[[tempdata]],[[data2]]]
        
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
      chk1=file_search(outdir+tar_name, count=nchk1)      
      if nchk1 eq 0 then spawn, 'mkdir -p '+outdir+tar_name

      ; check for other files with similar name (e.g. the ones with same size raster but extra row)
      ; because we don't want to overwrite these!
      ; check output directory exists and make if not
      check_again:
      fileout=outdir+tar_name+'/'+tar_name+'_'+xdim+'_'+ydim+'_p1_'+ind+'.sav'
      
      chk2=file_search(fileout, count=nchk2)
      if nchk2 gt 0 then begin
        ind2=ind+1
        ind=strtrim(ind2,2)
        ; check again
        goto, check_again
      endif 
      
      save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
        simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
        
    endfor  ; end loop over file
    
  endfor  ; end loop over target directory
  
  
  
  print, 'End of program'
  print, 'Concatenate files with concat_files.pro'
  
end