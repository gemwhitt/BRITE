pro fits_to_sav2

; PURPOSE: Read in .fits files, get header info + data (images)
;          Save data + info in .sav file - standard format
;          Save stats file
;
; OUTPUT: HDXXXXX.sav, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
;
;FOR PLOTTING
!p.background=cgcolor('white')
  
Compile_opt idl2  ; necessary to make integers a long (32 bit) value, rather than a short (int / 16 bit value)
  
sat='BA'      ; choose satellite
field='CENTAURUS' ; choose field
  
skipfiles=0  ;   ; number of fits files to skip over - if faulty
  
indir='~/BRITE/'+sat+'/'+field+'/data/raw_fits/rasters/'  ; location of fits files
outdir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/' ; location where save files will go
  
newtargets,indir,targets    ; call newtargets to get all HD names contained in this set of fits files
  
target=targets              ; string arrays with target names
nroi=n_elements(target)     ; number of unique ROIs to "unpack"
  
fitsfiles=file_search(indir+'*.fits', count=nfits)  ; locate and count ALL fits files
  
counter=0+skipfiles ; start a counter to move through each fits file
  
for roi=1, nroi-1 do begin
  
  print, 'Starting to build arrays for '+target[roi]
  
  statsfile=outdir+'stats/'+strcompress(target[roi], /remove_all)+'_roi.txt'
    
  repeat begin
    
    ; DATA
    data=mrdfits(fitsfiles[counter], target[roi], hdr, status=status, /silent)  ; get data for this ROI
    if status ne 0 then goto, next_file      ; -1=fail (i.e. no data)
    ;if status ne 0 then stop      ; -1=fail (i.e. no data)
    
    ; HEADER INFO
    ext0=mrdfits(fitsfiles[counter], 0, header, /silent) ; get the julian date
    
    ; TELEMETRY
    telem=mrdfits(fitsfiles[counter], 'TELEMETRY', telhead, /silent, status=stat1)
    if stat1 ne 0 then stop ; 0=successful completion
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    
    ; check name  
    extname=sxpar(hdr, 'EXTNAME')
    if extname ne target[roi] then stop         ; extname does not match target name
    
    ;check exposure time
    exposure_ttl=(sxpar(header, 'EXP_TTL'))/1000.         ; check total exposure time
    if exposure_ttl ne 1 then GOTO, next_file             ; if file is stacked then skip over - FOR NOW
      
    extname=strcompress(extname, /remove_all) ; take spaces out extname
      
    if n_elements(data1) eq 0 then begin  ; first entry
      
      data1=data
      
      sdata1=size(data1, /dim)
      
      jd=sxpar(header, 'JD-OBS')
      exp_num=sxpar(header, 'EXP_NUM')
      exp_time=(sxpar(header, 'EXP_TIME'))/1000.
      exp_ttl=exposure_ttl    ; in seconds
      
      ccd_temp=[telem.temp0,telem.temp1,telem.temp2,telem.temp3]
      stemp=strtrim(fix(average(ccd_temp)),2)
      
      roi_loc=[sxpar(hdr, 'ROI_X1'),sxpar(hdr, 'ROI_X2'),sxpar(hdr, 'ROI_Y1'),sxpar(hdr, 'ROI_Y2')]
      ra_dec=[sxpar(hdr, 'NOM_RA'),sxpar(hdr, 'NOM_DEC')]
      roi_name=[extname]
      
      caldat, jd, mon, day, yr
      sdate=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
      
      ; save stats
      openw,lun, statsfile, /get_lun, /append
      printf,lun, extname, sat, field, 'start_date', sdate, 'start_temp', stemp, $
        'raster_size', strtrim(sdata1[0],2), 'by', strtrim(sdata1[1],2), 'at', strtrim(counter,2), $
          format='(a10,x,a3,x,a10,x,a10,x,a10,x,a10,x,a3,x,a12,x,a3,x,a3,x,a3,x,a3,x,a10)'
      free_lun, lun
      
    endif else begin
      
      ; check dimensions agree
      snew=size(data, /dim)       ; size of new array
        
      sdata1=size(data1, /dim)      ; size of previous arrays
        
      if snew[0] eq sdata1[0] AND snew[1] eq sdata1[1] then begin ; append arrays
        
        data1=[[[data1]],[[data]]] 
        
        jd=[jd,sxpar(header, 'JD-OBS')]
        exp_num=[exp_num,sxpar(header, 'EXP_NUM')]
        exp_time=[exp_time,(sxpar(header, 'EXP_TIME'))/1000.]
        exp_ttl=[exp_ttl,exposure_ttl]    ; in seconds
        
        ccd_temp=[[ccd_temp],[telem.temp0,telem.temp1,telem.temp2,telem.temp3]]
        
        roi_loc=[[roi_loc],[sxpar(hdr, 'ROI_X1'),sxpar(hdr, 'ROI_X2'),sxpar(hdr, 'ROI_Y1'),sxpar(hdr, 'ROI_Y2')]]
        ra_dec=[[ra_dec],[sxpar(hdr, 'NOM_RA'),sxpar(hdr, 'NOM_DEC')]]
        roi_name=[roi_name,extname]
        
        if skipfiles gt 0 then begin
          openw,lun, statsfile, /get_lun, /append
          printf, lun, 'skipped', strtrim(skipfiles,2), 'files', format='(a8,x,a10,x,a6)'
          free_lun, lun
          
          skipfiles=0
        endif
        
      endif else begin
        
        skipfiles=skipfiles+1
        
        if skipfiles gt 1 then goto, next_file
        
        ; record stats
        openw,lun, statsfile, /get_lun, /append
        printf, lun, 'change_from', strtrim(sdata1[0],2), 'by', strtrim(sdata1[1],2), 'to', $
            strtrim(snew[0],2), 'by', strtrim(snew[1],2), 'at', strtrim(counter,2), $
            format='(a12,x,a3,x,a3,x,a3,x,a3,x,a3,x,a3,x,a3,x,a3,x,a10)'
          free_lun, lun
            
        goto, next_file
        
      endelse
        
    endelse

    next_file:
      
    counter=counter+1      ; increase counter by 1
      
  endrep until counter eq nfits  ; repeat until all ROIs have been processed
    
  ;check for duplicate times
  utime=jd[uniq(jd, sort(jd))]
  if n_elements(utime) lt n_elements(jd) then stop
  
  nsaved=strtrim(n_elements(jd),2)
  
  ; write out end stats?
  openw,lun, statsfile, /get_lun, /append
  printf, lun, nsaved, 'files_saved_out_of', strtrim(nfits,2), format='(a10,x,a20,x,a10)'
  free_lun, lun
    
  fileout=outdir+extname+'_p0.sav' ; make file
    
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl
    
  print, 'Completed .sav file for '+target[roi]
    
  ;reset counter
  undefine, jd, roi_name, exp_num, ra_dec, data1, roi_loc, ccd_temp, exp_time, exp_ttl
  skipfiles=0
  
  counter=0+skipfiles
  
    
endfor  ; end loop over this roi
  
print, 'end time, ',systime()
print, 'End of program'
print, 'Run SIMBAD info '+outdir
  
end