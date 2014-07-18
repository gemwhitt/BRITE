pro fits_to_sav3b

  ; PURPOSE: Read in .fits files, get header info + data (images)
  ;          Save data + info in .sav file - standard format
  ;          
  ; Modified from fits_to_sav3.pro - only save up to 500 images at a time - use target directories and an additional counter
  ;
  ; Run AFTER roi_fitsinfo.pro - use the info files to read only the necessary fits files - ignore the others
  ;
  ; OUTPUT: HDXXXXX/HDXXXXX_j.sav, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
  ;         where j increments after 500 images have been saved
  ;
  Compile_opt idl2  ; necessary to make integers a long (32 bit) value, rather than a short (int / 16 bit value)
  
  sat='BA'      ; choose satellite
  field='ORION' ; choose field
  
  fitsdir='~/BRITE/'+sat+'/'+field+'/data/raw_fits/rasters/'  ; location of fits files
  savdir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'            ; location where save files will go
  txtdir='~/BRITE/'+sat+'/'+field+'/data/raw_txt/'            ; location of txt info files
  
  target_files=file_search(txtdir+'HD*.txt', count=nroi)  ; locate the target (info) files
  targets=file_basename(target_files, '_fitsinfo.txt')    ; get the target names
  
  for roi=0, nroi-1 do begin                              ; loop over number of targets/rois
  
    ; create folder for this ROI - if one already exists then stop and check
    roi_dir=savdir+targets[roi]
    chk=file_search(roi_dir, count=nchk)
    if nchk eq 0 then spawn, 'mkdir '+roi_dir
    
    ; read in the info file and determine number of unique parameters - based on the last 4 variables
    readcol, target_files[roi], fname, exp_time, exp_ttl, xdim, ydim, $
      format='(x,x,x,a,x,x,f,f,i,i)'
    
    exp_time=fix(exp_time*10.)
    exp_ttl=fix(exp_ttl*10.)
    
    ; combine the observation params into one string
    obs_param=strtrim(exp_time,2)+'_'+strtrim(exp_ttl,2)+'_'+strtrim(xdim,2)+'_'+strtrim(ydim,2)
    
    undefine, xdim, ydim
    
    uparam=obs_param[uniq(obs_param, sort(obs_param))]    ; get unique parameter combination
    nuniq=n_elements(uparam)                              ; this is how many output files will be made
    ; ... for this target
    
    for obs=0, nuniq-1 do begin                          ; loop over each unique parameter combination
    
      ; get indices in obs_param to locate fits files which belong to this group of observation
      ii=where(obs_param eq uparam[obs], nfits)
      
      counter=0     ; this is to count through nfits - for these obsevration parameters
      j=0           ; this is to count 5000 files and then reset
      saved=0       ; this index will be denoted in the filename
      
      print, 'Starting to process '+targets[roi]+' at '+strtrim(systime(),2)
      
      repeat begin  ; loop over fits files
      
        fitsfile=fitsdir+fname[ii[counter]]+'.fits'
        
        ; get fitsfile info to get extnames
        fits_info, fitsfile, extname=extnames, /silent
        
        ; remove whitespace from extnames
        extnames1=strcompress(extnames, /remove_all)
        
        ; match targets[roi] to extnames1
        xx=where(extnames1 eq targets[roi], nmatch)
        
        if nmatch ne 1 then stop  ; if either there is no match - or > 1 match stop
        
        ; DATA
        data=mrdfits(fitsfile, extnames[xx], hdr, status=status, /silent)  ; get data for this ROI
        if status ne 0 then stop
        
        ; HEADER INFO
        ext0=mrdfits(fitsfile, 0, header, /silent) ; get the julian date etc
        undefine, ext0
        
        ; TELEMETRY
        telem=mrdfits(fitsfile, 'TELEMETRY', telhead, /silent, status=stat1)
        if stat1 ne 0 then stop ; 0=successful completion
        undefine, telhead
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        ; check name
        extname=sxpar(hdr, 'EXTNAME')
        if extname ne extnames[xx] then stop         ; extname does not match target name
        
        if n_elements(data1) eq 0 then begin  ; first entry
          data1=data
          
          jd=sxpar(header, 'JD-OBS')
          exp_num=sxpar(header, 'EXP_NUM')
          exp_time=(sxpar(header, 'EXP_TIME'))/1000.
          exp_ttl=(sxpar(header, 'EXP_TTL'))/1000.         ; total exposure time
          undefine, header
          
          ccd_temp=[telem.temp0,telem.temp1,telem.temp2,telem.temp3]
          undefine, telem
          
          roi_loc=[sxpar(hdr, 'ROI_X1'),sxpar(hdr, 'ROI_X2'),sxpar(hdr, 'ROI_Y1'),sxpar(hdr, 'ROI_Y2')]
          ra_dec=[sxpar(hdr, 'NOM_RA'),sxpar(hdr, 'NOM_DEC')]
          roi_name=extname
        endif else begin
          data1=[[[data1]],[[data]]]
          
          jd=[jd,sxpar(header, 'JD-OBS')]
          exp_num=[exp_num,sxpar(header, 'EXP_NUM')]
          exp_time=[exp_time,(sxpar(header, 'EXP_TIME'))/1000.]
          exp_ttl=[exp_ttl,(sxpar(header, 'EXP_TTL'))/1000.]    ; in seconds
          
          ccd_temp=[[ccd_temp],[telem.temp0,telem.temp1,telem.temp2,telem.temp3]]
          
          roi_loc=[[roi_loc],[sxpar(hdr, 'ROI_X1'),sxpar(hdr, 'ROI_X2'),sxpar(hdr, 'ROI_Y1'),sxpar(hdr, 'ROI_Y2')]]
          ra_dec=[[ra_dec],[sxpar(hdr, 'NOM_RA'),sxpar(hdr, 'NOM_DEC')]]
          roi_name=[roi_name,extname]
        endelse
        
        undefine, data
        undefine, hdr
        
        counter=counter+1      ; increase counter by 1
        j=j+1
        
        if j eq 4999 then begin
          
          ;check for duplicate times
          utime=jd[uniq(jd, sort(jd))]
          if n_elements(utime) lt n_elements(jd) then stop
          
          fileout=roi_dir+'/'+targets[roi]+'_'+uparam[obs]+'_p0_'+strtrim(saved,2)+'.sav' ; make file
          save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl
          
          ;reset j
          j=0
          saved=saved+1
          undefine, jd, roi_name, exp_num, ra_dec, data1, roi_loc, ccd_temp, exp_time, exp_ttl
        endif
        
      endrep until (counter eq nfits)  ; repeat until all ROIs have been processed
      
      ;check for duplicate times
      utime=jd[uniq(jd, sort(jd))]
      if n_elements(utime) lt n_elements(jd) then stop
      
      fileout=roi_dir+'/'+targets[roi]+'_'+uparam[obs]+'_p0_'+strtrim(saved,2)+'.sav' ; make file
      save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl
      
      ;reset counter
      counter=0
      j=0
      saved=0
      undefine, jd, roi_name, exp_num, ra_dec, data1, roi_loc, ccd_temp, exp_time, exp_ttl
      
      print, 'Saved output file for '+targets[roi]+', unique file '+strtrim(obs+1,2)+' out of '+strtrim(nuniq,2)
      print, 'at ', systime()
      
    endfor  ; end loop over this unique observation sequence
    
  endfor ; end loop over this target
  
  spawn, 'sudo purge'
  
  print, 'end time, ',systime()
  print, 'End of program'
  print, 'Run SIMBAD info '+savdir
  
end