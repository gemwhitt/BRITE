pro roi_make_sav

; Program to read in all data from a single observation window and make a save file with the required info
; - if ROI has already been observed then add to this save file
;
;; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')
;
Compile_opt idl2  ; necessary to make integers a long (32 bit) value, rather than a short (int / 16 bit value)
  
sat='UB'  ; swap betweem UB or BA
  
file_id='*/' ; some string specifying which file or group of files to include
  
indir=file_search('~/BRITE/data/'+sat+'/fits_data/rasters/ORION_1013to0314/'+file_id, count=ndir)  ; location of fits files
  
outdir='~/BRITE/data/'+sat+'/roi_raw_sav/ORION/' ; location where save files will go
  
; loop over each observation window, make .save files for each roi in each obs window
for i=1, ndir-1 do begin
  
  print, 'start time  ', systime()
  
  counter=0
  
  dirname=strmid(file_basename(indir[i]),0) ; extract name to append to output files
  print, dirname
  
  ; establish number of unique ROI's in these fits files
  fitsfiles=file_search(indir[i]+'/*.fits', count=nfits)
  
  next_part:
  data=mrdfits(fitsfiles[counter], 0, header, /silent)
  
  nroi0=sxpar(header, 'NUM_ROI') 
  
  ; cross-check number of rois in all fits files match the first
  ; if part b:
  
  for j=counter, nfits-1 do begin
    
    data=mrdfits(fitsfiles[j], 0, header, /silent)
    nroi=sxpar(header, 'NUM_ROI') 
    
    if nroi ne nroi0 then begin ; some files have different number of extensions - split the files accordingly
      
      print, 'Some extensions do not match'
      print, 'Splitting the data at '+strtrim(j-1,2)
      
      split=j
      stop
      goto, start
         
    endif else split=nfits
    
  endfor
  
  start:
  for k=1, nroi0 do begin  ; i tracks each roi (target)
      
    repeat begin
      
      ; first get the julian date from the 0th extension header
      data0=mrdfits(fitsfiles[counter], 0, header, /silent)
      jul_date=sxpar(header, 'JD-OBS')
      exposure_number=sxpar(header, 'EXP_NUM')
      exposure_time=(sxpar(header, 'EXP_TIME'))/1000.  ; in seconds
      exposure_ttl=(sxpar(header, 'EXP_TTL'))/1000.    ; in seconds
           
      ; then get the temperature values from the telemetry binary table
      telem=mrdfits(fitsfiles[counter], nroi0+2, header, /silent, status=stat1)
          
      if stat1 ne 0 then avg_temp='Nan' else begin
        
      temp0=telem.temp0
      temp1=telem.temp1
      temp2=telem.temp2
      temp3=telem.temp3
      
      ntemp=n_elements([temp0,temp1,temp2,temp3])
      
      avg_temp=(total([temp0,temp1,temp2,temp3]))/float(ntemp)
      
      endelse
   
      data=mrdfits(fitsfiles[counter], k, hdr, status=status, /silent)
      
      if status eq -1 then goto, next_file  ; no data for this ROI in this file
      
      ; i is the number of HDU's to skip from the current position - start with 1
      ; status is an integer number, 0=success, 2=end of file, -1=fail
      
      ; get some header info for this ROI
      naxis1=sxpar(hdr, 'NAXIS1')
      naxis2=sxpar(hdr, 'NAXIS2')
      extname=sxpar(hdr, 'EXTNAME')
      ; remove spaces from extname
      extname=strcompress(extname, /remove_all)
      
      ; if naxis eq 32 not 33 then add a dimension
      temp_data=lonarr(33, 33)
      temp_data[0:naxis1-1,0:naxis2-1]=data    
      data=temp_data
      
      x1=sxpar(hdr, 'ROI_X1')
      x2=sxpar(hdr, 'ROI_X2')
      y1=sxpar(hdr, 'ROI_Y1')
      y2=sxpar(hdr, 'ROI_Y2')
      ra=sxpar(hdr, 'NOM_RA')
      dec=sxpar(hdr, 'NOM_DEC')
;      stop
      ; calculate center of FOV
      xc=((x2-x1)/2.)+x1
      yc=((y2-y1)/2.)+y1
      
      ; build arrays
      if n_elements(jd) eq 0 then begin
        roi_name=extname
        exp_number=exposure_number
        roi_dim=[x1,x2,y1,y2]
        ra_dec1=[ra,dec]
        jd=jul_date
        data1=data
        ccd_temp=avg_temp
        exp_time=exposure_time
        exp_ttl=exposure_ttl
            
      endif else begin
        
        ; check this is the same extension
        if extname ne roi_name[0] then stop
        
        roi_name=[roi_name,extname]
        exp_number=[exp_number, exposure_number]
        ra_dec1=[[ra_dec1], [ra,dec]]
        jd=[jd, jul_date]
        data1=[[[data1]], [[data]]]
        ccd_temp=[ccd_temp,avg_temp]
        exp_time=[exp_time,exposure_time]
        exp_ttl=[exp_ttl,exposure_ttl]
        
      endelse
      
      next_file:
      ; increase counter by 1 
      counter=counter+1
      
      ; repeat until all ROIs have been processed
    endrep until counter eq split
    
    if split lt nfits then goto, next_part  ; repeat for files with different number ROI
    
    ; make file
    fileout=outdir+dirname+'_'+strtrim(extname,2)+'_p0.sav'
        
      exp_num=exp_number
      
      save, filename=fileout, roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl
      
      ;reset counter
      counter=0
      undefine, jd

  endfor  ; end loop over this roi
;  stop
endfor ; end loop over this directory (observation window)
  
  print, 'end time',systime()
  print, 'End of program'
  print, 'see output files in '+outdir
  
end