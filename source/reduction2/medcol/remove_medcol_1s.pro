pro remove_medcol_1s

; Do median column removal of raw_sav files using the full frame (medcol1) or raster (medcol3) values
; - decided upon by analyse_medcols.pro - restore output file
;
; Outout: save new file in p2 - with medcols removed
; 
Compile_opt idl2
  
sat='UB'
field='ORION'
   
indir='~/BRITE/'+sat+'/'+field+'/data/p1/'; CONTAINS ALL RASTER DIRECTORIES

outdir='~/BRITE/'+sat+'/'+field+'/data/p2/'; CONTAINS ALL RASTER DIRECTORIES
    
results_file='~/BRITE/'+sat+'/'+field+'/reduction/medcols/medcol_results.txt'
readcol, results_file, fname, medres1, medres3, sigres1, sigres3, result, format='a,d,d,d,d,a)', skipline=1

nf=n_elements(fname)

for ff=0, nf-1 do begin  ; begin loop over each file
  
  filein=file_search(indir+'HD*/'+fname[ff]+'.sav')
  
  filename=file_basename(filein)
   
      ; restore the whole file - because we are re-saving
      restore, filein   ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol1, medcol2, ndead, nnlin, medcol3, flag, pix50
         
      nimg=n_elements(jd)
      
      tempdata=data1*0
      
      s=size(data1, /dim)
      
      ; which medcol to use?
      if result[ff] eq 'medcol1' then medcol=medcol1
      if result[ff] eq 'medcol3' then medcol=medcol3
      
      ; loop over each image and each column -
      ;   subtract the medcol1/medcol2 value, make any negative pixels non-zero
      for im=0, nimg-1 do begin
        
        im0=data1[*,*,im]
        
;        wset, 0
;        plot_image, bytscl(im0, 0, 500)
                     
        for c=0, s[0]-1 do im0[c,*]=im0[c,*]-medcol[c,im]
                
        ; Check for negative pixels
        neg=where(im0 lt 0, nneg)
        
        if nneg gt 0 then begin
        
          ineg=array_indices(im0, neg)
          
           ;replace negative pixels with a positive integer, not zero
          for nn=0, nneg-1 do im0[ineg[0,nn],ineg[1,nn]]=0
          
        endif
        
;        wset, 1
;        plot_image, bytscl(im0, 0, 500)
;        stop
        
        ; save new data2 in tempdata
        tempdata[*,*,im]=im0
      
      endfor  ; end loop over image - keep going until all images are processed
      
      ; rename tempdata to data1 and save result in new file
      data1=tempdata
      
      fileout=outdir+filename
          
      save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
        
    endfor  ; end loop over file
      
  
  
  print, 'End of program'
  print, 'Concatenate files with concat_files.pro'
  
end