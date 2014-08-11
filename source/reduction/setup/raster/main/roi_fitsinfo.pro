pro roi_fitsinfo

  ; PURPOSE: Read in .fits files, get header info + data (images)
  ;          Save data + info in .txt file for each target - append this for each fits file
  ;          Convert .txt files to .sav files using another program
  ;          Save stats file
  ;
  ; Modified from fits_to_sav2.pro on 16th July
  ;
  ; OUTPUT: HDXXXXX.txt, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
  ;
  ;FOR PLOTTING
  !p.background=cgcolor('white')
  
  Compile_opt idl2  ; necessary to make integers a long (32 bit) value, rather than a short (int / 16 bit value)
  
  sat='BA'      ; choose satellite
  field='CENTAURUS' ; choose field
  
  skipfiles=0  ;   ; number of fits files to skip over - if faulty
  
  dates=['2014-06','2014-07']
  
  indir='~/BRITE/'+sat+'/'+field+'/data/raw_fits/rasters/'+dates+'/'  ; location of fits files
  outdir='~/BRITE/'+sat+'/'+field+'/data/raw_txt/2014-0607/' ; location where save files will go
  
  fitsfiles=file_search(indir+'*.fits', count=nfits)  ; locate and count ALL fits files
  
  badfiles=outdir+'status.txt' ; create badfiles.txt to note any .fits files which have no data or are corrupt
  
  print, 'Making fits_info files for '+sat+' '+field
  print, 'Number of fits files: '+strtrim(nfits,2)
  print, 'Start time is: '+systime()
  
  for f=skipfiles, nfits-1 do begin
    
    ; calculate percentage of fits done
    pcent=float(f)/float(nfits)*100.
  
    fname=file_basename(fitsfiles[f], '.fits')
    
    if f eq 0. then print, 'Processing '+fname+' 0th file'
    if f mod 100. eq 0 then print, strtrim(pcent,2)+' percent done'
    
    ; check that file is not corrupt, else skipover
    ; HEADER INFO
    ext0=mrdfits(fitsfiles[f], 0, header, /silent, status=status) ; get the julian date
    
    if status ne 0 then begin
      fmt1='(a5,x,a10,x,a'+strtrim(strlen(fname),2)+',x,a12)'
      
      stop
      openw, lun, badfiles, /get_lun, /append
      printf, lun, sat, field, fname, 'empty_file', format=fmt1
      free_lun, lun
      
      goto, next_file
    endif
    
    ; first check that .fits file has target data - if not print this in status file
    fits_info, fitsfiles[f], /silent, extname=extnames
    
    rois=where(strmid(extnames,0,2) eq 'HD', n_roi)
    
    if n_roi le 0 then begin
      fmt1='(a5,x,a10,x,a'+strtrim(strlen(fname),2)+',x,a12)'
      
      openw, lun, badfiles, /get_lun, /append
      printf, lun, sat, field, fname, 'no_targets', format=fmt1
      free_lun, lun
      
      goto, next_file
    endif
    
    ; check that telemetry info exists
    ; TELEMETRY
    telem=mrdfits(fitsfiles[f], 'TELEMETRY', telhead, /silent, status=stat1)
    if stat1 ne 0 then begin
      fmt1='(a5,x,a10,x,a'+strtrim(strlen(fname),2)+',x,a12)'
      
      openw, lun, badfiles, /get_lun, /append
      printf, lun, sat, field, fname, 'no_telemetry', format=fmt1
      free_lun, lun
      
      goto, next_file
    endif
    
    
    ; if the program gets to this point, then the fits file is good and there are targets to track
    ; for each target with data in the fitsfile - update that target's .txt file with info on raster size and jd
    
    ; HEADER info
    jd=sxpar(header, 'JD-OBS')
    
    caldat, jd, mon, day, yr
    obsdate=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
    
    exp_num=sxpar(header, 'EXP_NUM')
    exp_time=(sxpar(header, 'EXP_TIME'))/1000.       ; in seconds
    exp_ttl=(sxpar(header, 'EXP_TTL'))/1000.         ; in seconds - if stacked
    
    ; check extention by extention and append target.txt files
    for roi=1, n_roi do begin
      data=mrdfits(fitsfiles[f], roi, hdr, status=status, /silent)  ; check header for this roi
      if status ne 0 then stop      ; -1=fail (i.e. no data)
      
      ; check name
      extname=sxpar(hdr, 'EXTNAME')
      if extname ne extnames[roi] then stop         ; extname does not match target name
      
      roi_name=strcompress(extname, /remove_all) ; take spaces out extname
      
      ; check dimensions of data array
      sdata=size(data, /dim)       ; size of new array
      
      ; write info to .txt file for this target
      fileout=outdir+roi_name+'_fitsinfo.txt'
      openw, lun, fileout, /get_lun, /append
      printf, lun, sat, field, obsdate, fname, roi_name, exp_num, exp_time, exp_ttl, sdata[0], sdata[1], $
        format='(a5,x,a10,x,a10,x,a40,x,a10,x,i6,x,f6.3,x,f6.3,x,i3,x,i3)'
      free_lun, lun
    endfor  ; end loop over this ROI
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    next_file:
    
  endfor  ; end loop over this fitsfile
  
  print, 'End time is: '+systime()
  print, 'End of program'
  print, 'Now extract roi data with fits_to_sav3b'
  
end