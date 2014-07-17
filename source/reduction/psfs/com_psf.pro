pro com_psf

; Purpose: Use region_psf results to calulate the COM for each target PSF and save info in the .sav file

Compile_opt idl2

sat='BA'
field='ORION'

indir='~/BRITE/data/'+sat+'/p2/'+field+'/class/sav/'

statsdir='~/BRITE/data/'+sat+'/p2/'+field+'/class/stats/'
statsfile=statsdir+'region_psf_stats.txt'

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin

  hdname=file_basename(filesin[f], '_p2.sav')
  
  restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
  ;medcol1, medimg0, ndead, nsat, simbad_radec, vmag, bmag, $
  ;parlax, otype, sptype, iflag, ccd_temp2, psf_loc, flx_reg
  
  nimg=n_elements(jd)
  
  jd1=jd-jd[0]
  
  xy_com=fltarr(2,nimg)
  
  for im=0, nimg-1 do begin ; begin loop over images
    
    if iflag[im] ne 2 then continue
  
    im0=reform(data1[*,*,im])
    
    s0=size(im0, /dim)
    
    ; Use psf_loc to get cutout
    x1=(psf_loc[0,im]-1 > 0)
    x2=(psf_loc[1,im]+1 < s0[0]-1)
    y1=(psf_loc[2,im]-1 > 0)
    y2=(psf_loc[3,im]+1 < s0[1]-1)
    
    cutout=im0[x1:x2,y1:y2]
    
    s=size(cutout, /dim) ; length of array in x and y
    
    imtot=total(cutout)
    
    x_cen=total(total(cutout, 2)*indgen(s[0]))/imtot
    y_cen=total(total(cutout, 1)*indgen(s[1]))/imtot
    
    xy_com[0,im]=x_cen
    xy_com[1,im]=y_cen
    
  endfor  ; end loop over images 
  
  ; re-save .save file with xy_com added
  
  save, filename=filesin[f], roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
  medcol1, medimg0, ndead, nsat, simbad_radec, vmag, bmag, $
  parlax, otype, sptype, iflag, ccd_temp2, psf_loc, flx_reg, xy_com
  
endfor ; end loop over file
 
print, 'end of program'
end