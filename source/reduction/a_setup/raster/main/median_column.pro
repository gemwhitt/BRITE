pro median_column

; Program to remove warm columns from the rasters - values for column medians are saved in the data - at the end of columns
; 
Compile_opt idl2
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

sat='BA'
field='ORION'

indir='~/BRITE/data/'+sat+'/roi_raw_sav/ORION/'

outdir='~/BRITE/data/UB/p1/ORION/'  

filesin=file_search(indir+'*.sav', count=nfiles)

for i=0, nfiles-1 do begin
  
  fname=file_basename(filesin[i], '_p0.sav')
  ;print, fname
  
  restore, filesin[i]   ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
                        ;simbad_radec, vmag, bmag, parlax, otype, sptype
                        
                
  nimg=(size(data1, /dim))[2]
  
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  medcols=fltarr(xdim-1, nimg)
  
  medimg0=fltarr(nimg)
  medimg1=fltarr(nimg)
  ndead=intarr(nimg)  ; number of dead pixels on image
  nsat=intarr(nimg)   ; number of pixels above a given threshold on the image
  
  tempdata=lonarr(xdim-1, ydim-1,nimg)
  
 ; window,0, xsize=600, ysize=550, xpos=1500, ypos=100
 ; for ll=0, 42 do begin
   ; plot_image, bytscl(data1[*,*,0], 20, 500)
   ; stop
 ;   wait, 0.5
 ; endfor
 ; 
;stop
  for j=0, nimg-1 do begin
    
    data2=data1[0:xdim-1,0:ydim-2,j]
   
    ; check upper left corner for a non-zero value
    chkpix=data1[0,ydim-1,j]

    ;if chkpix gt 0 then medcols[*,j]=data1[0:xdim-2,ydim-1,j] else begin
    
     ; calculate the median of each column 
     for k=0, xdim-2 do medcols[k,j]=median(data2[k,0:ydim-2])
      
    ;endelse
   
    ; calculate image median
    medimg0[j]=median(data2)
    
    ; check for dead pixels
    dead=where(data2 le 0, num_dead)
    ndead[j]=num_dead
    
    ; check for saturated pixels
    sat=where(data2 ge 15000, num_sat)
    nsat[j]=num_sat
    
;    if num_sat gt 0 then stop
    
    ; subtract the median columns    
    for k=0, xdim-2 do tempdata[k,*,j]=data2[k,*]-medcols[k,j]
    
;    if num_sat gt 0 then stop
    
    ; make any negative values zero values instead
    data2=reform(tempdata[*,*,j])
    
    xx=where(data2 lt 0, nneg)
    
    if nneg gt 0 then begin
    
      loc2d=array_indices(data2, xx)
      
      data2[loc2d[0,*], loc2d[1,*]]=0
      
      xx=where(data2 lt 0, nneg)
      
      if nneg gt 0 then stop
      
      tempdata[*,*,j]=data2
      
    endif
    
    medimg1[j]=median(data2)
    
  endfor
  
  data1=tempdata
  
;  window,1, xsize=600, ysize=550, xpos=2300, ypos=100
;  plot_image, bytscl(data1[*,*,0], 20, 500)

xx=where(nsat eq 1024, nxx)
openw, lun, '~/Desktop/orion_sat_imgs.dat', /get_lun, /append
printf, lun, fname, nxx, format='(a15, x, i4)'
free_lun, lun

;cgminmax, data1

  ; save data
  ;fileout=outdir+fname+'_p1.sav'
  ;save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
  ;              medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype
               
;stop              
endfor

spawn, 'open '+ '~/Desktop/orion_sat_imgs.dat'


print, 'end of program'
end


