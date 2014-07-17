pro fits_to_sav

; PURPOSE: Read .fits files in a given folder, get header info + data (images) 
;          Save data + info in .sav file - standard format
;          
; OUTPUT: HDXXXXX.sav, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
;
;FOR PLOTTING
!p.background=cgcolor('white')

Compile_opt idl2  ; necessary to make integers a long (32 bit) value, rather than a short (int / 16 bit value)
  
sat='BA'      ; choose satellite
field='ORION' ; choose field

skipfiles=0  ;   ; number of fits files to skip over - if faulty
    
indir='~/BRITE/data/'+sat+'/fits_data/rasters/ORION/'  ; location of fits files  
outdir='~/BRITE/data/'+sat+'/roi_raw/ORION/' ; location where save files will go
  
newtargets,indir,targets    ; call newtargets to get all HD names contained in this set of fits files
  
target=targets              ; string arrays with target names
nroi=n_elements(target)     ; number of unique ROIs to "unpack" 

fitsfiles=file_search(indir+'*.fits', count=nfits)  ; locate and count ALL fits files

counter=0+skipfiles ; start a counter to move through each fits file
        
for roi=0, nroi-1 do begin  
    
  print, 'Starting to build arrays for '+target[roi]
  
  repeat begin
    
    data=mrdfits(fitsfiles[counter], target[roi], hdr, status=status, /silent)  ; get data for this ROI
    
    if status eq -1 then goto, next_file      ; -1=fail (i.e. no data)
    
    extname=sxpar(hdr, 'EXTNAME')
    
    if extname ne target[roi] then stop         ; extname does not match target name
    
    extname=strcompress(extname, /remove_all) ; take spaces out extname
    
    if n_elements(data1) eq 0 then data1=data else begin
      
      ; check dimensions agree
      sn=size(data, /dim)       ; size of new array
      
      sp=size(data1, /dim)      ; size of previous arrays
      
      if sn[0] eq sp[0] AND sn[1] eq sp[1] then data1=[[[data1]],[[data]]] else begin
        
        ; save stats
        statsfile=outdir+'stats/'+extname+'_roi.txt'
        openw,lun, statsfile, /get_lun, /append
        printf,lun, extname, sat, field, 'raster_size', strtrim(sp[0],2), 'by', strtrim(sp[1],2), 'change_to', $
          strtrim(sn[0],2), 'by', strtrim(sn[1],2), $
          'at', strtrim(counter,2), format='(a10,x,a3,x,a10,x,a12,x,a3,x,a3,x,a3,x,a10,x,a3,x,a3,x,a3,x,a3,x,a10)'
        free_lun, lun
        
        goto, next_file
      endelse
      
    endelse
  
      
  data0=mrdfits(fitsfiles[counter], 0, header, /silent) ; get the julian date 
        
  exposure_ttl=(sxpar(header, 'EXP_TTL'))/1000.         ; check total exposure time
  if exposure_ttl ne 1 then GOTO, next_file             ; if file is stacked then skip over - FOR NOW
        
  if n_elements(jd) eq 0 then begin                     ; start building arrays for these parameters
    jd=sxpar(header, 'JD-OBS')
    exp_num=sxpar(header, 'EXP_NUM')
    exp_time=(sxpar(header, 'EXP_TIME'))/1000.
    exp_ttl=exposure_ttl    ; in seconds
  endif else begin                                      ; append arrays for these parameters
    jd=[jd,sxpar(header, 'JD-OBS')]
    exp_num=[exp_num,sxpar(header, 'EXP_NUM')]
    exp_time=[exp_time,(sxpar(header, 'EXP_TIME'))/1000.]  
    exp_ttl=[exp_ttl,exposure_ttl]    ; in seconds
  endelse
        
  ; get the temperature values from the telemetry binary table
  telem=mrdfits(fitsfiles[counter], 'TELEMETRY', header, /silent, status=stat1)
        
  if stat1 ne 0 then stop ; 0=successful completion
        
  if n_elements(ccd_temp) eq 0 then $
    ccd_temp=[telem.temp0,telem.temp1,telem.temp2,telem.temp3] else $
    ccd_temp=[[ccd_temp],[telem.temp0,telem.temp1,telem.temp2,telem.temp3]]       
            
  ; get roi header info
  if n_elements(roi_loc) eq 0 then begin
    roi_loc=[sxpar(hdr, 'ROI_X1'),sxpar(hdr, 'ROI_X2'),sxpar(hdr, 'ROI_Y1'),sxpar(hdr, 'ROI_Y2')]
    ra_dec=[sxpar(hdr, 'NOM_RA'),sxpar(hdr, 'NOM_DEC')]
    roi_name=[extname]
  endif else begin
    roi_loc=[[roi_loc],[sxpar(hdr, 'ROI_X1'),sxpar(hdr, 'ROI_X2'),sxpar(hdr, 'ROI_Y1'),sxpar(hdr, 'ROI_Y2')]]
    ra_dec=[[ra_dec],[sxpar(hdr, 'NOM_RA'),sxpar(hdr, 'NOM_DEC')]]
    roi_name=[roi_name,extname]
  endelse
       
  xdim=sxpar(hdr, 'NAXIS1') ; dimension of raster in x
  ydim=sxpar(hdr, 'NAXIS2') ; dimension of raster in y
                    
   next_file:
   
   counter=counter+1      ; increase counter by 1
        
 endrep until counter eq nfits  ; repeat until all ROIs have been processed
 
 ;check for duplicate times
 utime=jd[uniq(jd, sort(jd))]
 if n_elements(utime) lt n_elements(jd) then stop

 fileout=outdir+'sav/'+extname+'_p0.sav' ; make file
      
 save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl
      
 print, 'Completed .sav file for '+target[roi]
           
 ;reset counter
 counter=0+skipfiles
 undefine, jd, roi_name, exp_num, ra_dec, data1, roi_loc, ccd_temp, exp_time, exp_ttl
      
endfor  ; end loop over this roi
 
print, 'end time, ',systime()
print, 'End of program'
print, 'Run SIMBAD info '+outdir
  
end