pro region_psf

; PURPOSE: Use label_region.pro to 1st. identify groups of pixels with DN > thr, 2nd. identify PSF pixels
;                                  3rd. measure borders of PSf for each image
;                                  
; Output: 1. Write out stats
;         2. re-save .sav files with PSF x and y limits
;                                  
; Follow-up with get_cen.pro                                
; 
; Date created: 14th June 2014
; Modified 17th JUne 2014
; 
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
                      ;parlax, otype, sptype, iflag, ccd_temp2
  
  psf_loc=roi_loc*0
  
  nimg=n_elements(jd)
  
  jd1=jd-jd[0]
  
  npix=intarr(nimg)
  flx=dblarr(nimg)
  
  for im=0, nimg-1 do begin ; begin loop over images
    
    im0=reform(data1[*,*,im])
    
    s=size(im0, /dim) ; length of array in x and y
    im1=lonarr(s+2)   ; append array with extra row/column (for label_region.pro)
    im1[1,1]=im0      ; ""
        
    thr=30
        
    r1=label_region(im1 gt thr, /all_neighbors)  ; find groups of pixels above thr
    r1=r1[1:s[0],1:s[1]]                         ; trim off border
   
    ; Get population and members of each blob:
    h=HISTOGRAM(r1, REVERSE_INDICES=ri, locations=loc)
    
    ; change the order in order of decreasing h
    neworder=reverse(sort(h))
    h2=h[neworder]
    loc2=loc[neworder]
    
    if n_elements(h) lt 2 then continue else begin
      npix[im]=h2[1]    ; number of pixels in PSF
      ind=loc2[1]   ; label of those pixels in r1 which correspond too the PSF
    endelse
      
    ; get x and y indices of PSF pixels using ind and r1 - use to measure borders and inital flux
    psfs=where(r1 eq ind, npsfs)
    ; check npsfs = h2[1]
    if npsfs ne h2[1] then stop
    xi=(array_indices(r1, psfs))[0,*]
    yi=(array_indices(r1, psfs))[1,*]
    ;plot_image, bytscl(r1, 0, max(loc))
    ;oplot, xi, yi, color=cgcolor('purple'), psym=2
    
    ; get borders of PSF
    psf_loc[*,im]=[min(xi),max(xi),min(yi),max(yi)]
    
    ; measure total DN in PSF pixels
    flx[im]=total(im0[xi,yi])
    
  endfor  ; end loop over images
  
  ; calculate median and sigma of variables - ignore 0 values
  keep=where(flx ne 0, nkeep)
  
  med_npix=median(npix[keep])
  sig_npix=robust_sigma(npix[keep]/med_npix)
  med_flx=median(flx[keep])
  sig_flx=robust_sigma(flx[keep]/med_flx)
  
  ; write out stats
  openw, lun, statsfile, /get_lun, /append
  printf, lun, roi_name[0], med_npix, sig_npix, med_flx, sig_flx, format='(a7,x,i3,x,f7.3,x,d14.2,x,f7.3)'
  free_lun, lun
  
  ; re-save .sav file with psf_loc & flx
  savout=filesin[f]
  
  flx_reg=flx
  
  save, filename=savout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    medcol1, medimg0, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype, iflag, ccd_temp2, psf_loc, flx_reg
  
endfor ; end loop over files

print, 'End of Program'
end

