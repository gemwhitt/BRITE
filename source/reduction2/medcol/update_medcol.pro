pro update_medcol 

; Purpose: Perform a preliminary medcol removal of the images using medcol2 (raster medcols)
; Find regions of interest (PSF and HPs) - convert these to NaN and then recalculate medcol raster 
; 
; Output: 
; 
; Variables Modified: None
; 
; Variables added: medcol3, flag, pix50
; 
; Where in processing: After save_medcol, before analyse_medcols
;
Compile_opt idl2

sat='UB'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/p1/' ; CONTAINS ALL FOLDERS FOR THE TARGETS

outdir=indir

filein=file_search(indir+'HD*/*.sav', count=nf)


for ff=0, nf-1 do begin  ; begin loop over each file
  
    ; restore the whole file - because we are re-saving
    restore, filein[ff]   ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol1, medcol2, ndead, nnlin
                                                    
    nimg=n_elements(jd)
    
    medcol3=medcol2*0
    flag=intarr(nimg)+2       ; all images flagged as 2 (target) - unless found otherwise
    pix50=float(data1*0)+1    ; array modelled on data1 to indicate all pixels above threshold
        
    ; loop over each image and each column -
    ;   subtract the medcol1/medcol2 value, make any negative pixels non-zero
    for im=0, nimg-1 do begin
    
      im0=data1[*,*,im]
      
      s=size(im0, /dim)
      
;      plot_image, bytscl(im0, 0, 500)
;      stop
      
      ; check that the number of columns in medcol2 and im0 agree
      sm2=size(medcol2, /dim)
      
      if s[0] ne sm2[0] then stop
      
      im1=im0*0       ; define im1 to be im0 with medcol2 removed
      for c=0, s[0]-1 do im1[c,*]=im0[c,*]-medcol2[c,im]
      
;      plot_image, bytscl(im1, 0, 500)
;      stop

      ; now find regions of interest in im1 with label_region
      ; first, add an extra row and column to the image
      im2=lonarr(s[0]+2,s[1]+2)
      im2[1,1]=im1
      
;      wset, 0
;      plot_image, bytscl(im2, 0, 500)
      
      ; make all negative values 0
      neg=where(im2 lt 0, nneg)
      if nneg gt 0 then begin
        neg2d=array_indices(im2, neg)
        
        for i=0, nneg-1 do im2[neg2d[0,i],neg2d[1,i]]=0
        
;        wset, 1
;        plot_image, bytscl(im2, 0, 500)
;        stop
      endif
      
      rim2=label_region(im2 ge 50, /all_neighbors)
      rim1=rim2[1:s[0],1:s[1]]
      
      hist=histogram(rim1, binsize=1, locations=loc, reverse_indices=ri)
      nhist=n_elements(hist)
      
      if nhist eq 1 then begin $  ; no regions brighter than 50 DN - flag=0 and continue
        flag[im]=0
        goto, next_image
      endif
      
      ; make all but the first group of pixels NaN values
      bright_ind=ri[ri[1]:ri[nhist]-1]  ; this is the index locations of pix50
      
      ; check that the number of pixels adds up
      if n_elements(bright_ind)+hist[0] ne s[0]*s[1] then stop
      
      ; make all illuminated pixels in im0 NaN       
;      wset, 0
;      plot_image, bytscl(im0, 0, 500)
;      stop
      
      ind2d=array_indices(im0, bright_ind)
      
      for pix=0, n_elements(bright_ind)-1 do pix50[ind2d[0,pix], ind2d[1,pix],im]=!Values.F_NAN
      
;      ; check
;      plot_image, bytscl(im1, 0, 500)
;      oplot, ind2d[0,*], ind2d[1,*], color=cgcolor('green'), psym=2
;      stop
      
      for pix=0, n_elements(bright_ind)-1 do im0[ind2d[0,pix], ind2d[1,pix]]=!Values.F_NAN
      
;      wset, 1
;      plot_image, bytscl(im0, 0, 500)
;      stop
      
;      ; check
;      plot_image, bytscl(im1, 0, 500)
;      oplot, ind2d[0,*], ind2d[1,*], color=cgcolor('green'), psym=2
;      stop

      ; now compute new medcol values
      for c=0, s[0]-1 do medcol3[c,im]=median(im0[c,*])
     next_image:
    endfor  ; end of loop over images
    
    ; now save the original file with medcol3 appended
    save, filename=filein[ff], roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol1, medcol2, ndead, nnlin, medcol3, flag, pix50
    
  endfor  ; end loop over file

print, 'End of Program'
print, 'Now run analysis followed by remove_medcol_1s with choice of medcol value (i.e. 1, 2, or 3)'
end