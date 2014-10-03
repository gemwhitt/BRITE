pro save_medcol

;   Calculate and record raster median columns (e.g. 24x24), record full-frame median columns (~2700x24)

;   Save medcol1 (full-frame) and medcol2 (raster) AND 
;   ... ndead (number of 0 values pixels) and nnlin (number of pixels > thr=9000)
;   
;   Variables Modified: data1 - trimmed if there is an extra row with the CCD medcol values
;   
;   Variables Added: medimg0, ndead, nnlin, medcol1, medcol2
;   
Compile_opt idl2

sat='UB'
field='CENTAURUS'

indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/2014_0607/'

tardir=file_search(indir+'HD*/', count=ntar) ; target directories

outdir='~/BRITE/'+sat+'/'+field+'/data/p1/2014_0607/'+file_basename(tardir)


for tar=0, ntar-1 do begin
  
  filesin=file_search(tardir[tar]+'/*.sav', count=nf) 
  
  for ff=0, nf-1 do begin
  
    ; restore the whole file
    restore, filesin[ff]   ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype
    
    ; get hdname and index for output file
    fname=file_basename(filesin[ff], '.sav')
    dashpos=strsplit(fname, '_')
    hdname=strmid(fname, 0, dashpos[1]-1)
    index=strmid(fname, dashpos[6], 1)
    
    ; if it's not a 1s non-stacked exposure then skip over (for now)
    if exp_ttl[0] ne 1. OR exp_time[0] ne 1. then continue
    
    nimg=n_elements(jd)
    
    ; define array for new data1
    newdata1=[]
    
    s=size(data1, /dim)
    
    medcol1=lonarr(s[0], nimg)  ; array of full-frame median columns - calculated onboard
    medcol2=lonarr(s[0], nimg)  ; array of raster-calculated column medians
    
    medimg0=fltarr(nimg)        ; array of whole image medians (raw)
    
    ndead=intarr(nimg)  ; initial number of dead pixels on image
    nnlin=intarr(nimg)   ; number of pixels above a given threshold on the image
    
    for im=0, nimg-1 do begin
    
      ; first check that top row is medcol values and not just data
      ; rasters will usually contain an even number of rows/columns - ASSUMPTION!!!
      ; so check this modify accordingly
      
      if s[1] mod 2 eq 0 then data2=data1[0:s[0]-1,0:s[1]-1,im] $ ; even number of rows - i.e rasters don't contain extra info
        else begin                                                ; odd number of rows - i.e rasters contain medcol values
          data2=data1[0:s[0]-1,0:s[1]-2,im]    
          
          ; extract the fullframe median columns
          medcol1[*,im]=data1[*,s[1]-1,im]
        endelse
      
      ; calculate median from raster
      s2=size(data2, /dim)
      for k=0, s2[0]-1 do medcol2[k,im]=median(data2[k,*])
                     
      medimg0[im]=median(data2)           ; calculate image median
      
      dead=where(data2 le 0, num_dead)    ; check for dead pixels
      ndead[im]=num_dead
      
      nonlin=where(data2 ge 9000, num_nlin)  ; check for nonlinear
      nnlin[im]=num_nlin
      
      newdata1=[[[newdata1]],[[data2]]]
      
    endfor  ; end loop over this image
    
    ; resave newdata1 as data1
    data1=newdata1
    
    ; check output directory exists
    chk=file_search(outdir[tar],count=nchk)
    if nchk eq 0 then spawn, 'mkdir -p '+outdir[tar]
      
    ; save output
    fileout=outdir[tar]+'/'+hdname+'_'+strtrim(s2[0],2)+'_'+strtrim(s2[1],2)+'_'+index+'.sav'
  
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol1, medcol2, ndead, nnlin
   
  
endfor  ; end loop over this file
  
endfor  ; end loop over this target

print, 'End of program'
print, 'Do update_medcol'

end
