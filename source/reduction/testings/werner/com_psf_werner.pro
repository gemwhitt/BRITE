pro com_psf_werner

  ; Purpose: Calulate the COM for each target PSF and save info in the .sav file
  
  Compile_opt idl2
  
  indir='~/BRITE/TESTSETS/werner4lc/p2/subsets/'
   
  filesin=file_search(indir+'*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, nf-1 do begin
  
    hdname=file_basename(filesin[f], '_p2.sav')
    
    restore, filesin[f] ;jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag, flag, psf_loc
    
    nimg=n_elements(jd)
    
    jd1=jd-jd[0]
    
    xy_com=fltarr(2,nimg)
    
    for im=0, nimg-1 do begin ; begin loop over images
    
      if flag[im] ne 2 then continue
      
      im0=reform(data1[*,*,im])
      
      s0=size(im0, /dim)
      
      ; Use psf_loc to get cutout
;      x1=(psf_loc[0,im]-1 > 0)
;      x2=(psf_loc[1,im]+1 < s0[0]-1)
;      y1=(psf_loc[2,im]-1 > 0)
;      y2=(psf_loc[3,im]+1 < s0[1]-1)
      x1=(psf_loc[0,im])
      x2=(psf_loc[1,im])
      y1=(psf_loc[2,im])
      y2=(psf_loc[3,im])
      
      if (x2-x1)*(y2-y1) lt 50 then begin
        flag[im]=0
        goto, next_im
      endif
      
      cutout=im0[x1:x2,y1:y2]
      
      s=size(cutout, /dim) ; length of array in x and y
      
      imtot=total(cutout)
       
      x_cen=total(total(cutout, 2)*indgen(s[0]))/imtot
      y_cen=total(total(cutout, 1)*indgen(s[1]))/imtot
      
      xy_com[0,im]=x_cen+x1
      xy_com[1,im]=y_cen+y1
      
      ; plot results to check
      ;plot_image, bytscl(im0, 20, 200)
      ;oplot, x1+[xy_com[0,im]], y1+[xy_com[1,im]], color=cgcolor('orange'), psym=2
      ;stop
      next_im:
    endfor  ; end loop over images
    
    ; re-save .save file with xy_com added
    
    save, filename=filesin[f], jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag, flag, psf_loc, xy_com
      
  endfor ; end loop over file
  
  print, 'end of program'
end