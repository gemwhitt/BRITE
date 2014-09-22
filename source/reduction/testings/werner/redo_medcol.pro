pro redo_medcol, im0, psf, medcol2, medcol3, target, med_diff

; program for re-doing the initial median column removal - avoing target PSFs and HPs
;
Compile_opt idl2

thr=50

indir='~/BRITE/TESTSETS/werner4lc/p3/'
outdir='~/BRITE/TESTSETS/werner4lc/p4/'

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin

  fname=file_basename(filesin[f], '.sav')
  
  restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
  ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1, flag, $
  ;psf_loc, npix_psf, hp_loc, nhp, modelpsf
  
  newmodel=modelpsf*0
  
  nimg=n_elements(jd)
  
  jd1=jd-jd[0]
  
  xy_com1=fltarr(2,nimg)
  xy_com2=fltarr(2,nimg)
  medcol3=medcol2*0
  median_diff=intarr(nimg)
  
  for im=0, nimg-1 do begin ; begin loop over images
  
    if flag[im] ne 2 then continue
    
    im0=reform(data1[*,*,im])
    
    s0=size(im0, /dim)
    
    ; check if any HPs remaining and neutralise these...
    ihp=where(hp_loc[*,*,im] eq 1, nhp2)
    for hp=0, nhp[im]-1 do begin
    
      if nhp[im] ne nhp2 then stop
      
      ; make a cutout around each HP
      xhp=(array_indices(reform(hp_loc[*,*,im]), ihp))[0,hp]
      yhp=(array_indices(reform(hp_loc[*,*,im]), ihp))[1,hp]
      
      ;        wset, 0
      ;        plot_image, bytscl(im0, 0, 50)
      ;        oplot, [xhp], [yhp], color=cgcolor('purple'), psym=2
      
      ;        plot_image, bytscl(hp_loc[*,*,im], 0, 1)
      ;        oplot, [xhp], [yhp], color=cgcolor('purple'), psym=2
      
      ; make a subsquare around the HP providing that it's not at the edge of the image
      if xhp eq 0 then x2=xhp+2 else x2=(xhp+1 < (s0[0]-1))
      if xhp eq s0[0]-1 then x1=xhp-2 else x1=((xhp-1) > 0)
      if yhp eq 0 then y2=yhp+2 else y2=(yhp+1 < (s0[1]-1))
      if yhp eq s0[1]-1 then y1=yhp-2 else y1=((yhp-1) > 0)
      
      ;        wset, 1
      subsq=im0[x1:x2,y1:y2]
      ;        plot_image, bytscl(subsq, 0, 200)
      
      im0[xhp,yhp]=median(subsq)
      data1[xhp,yhp,im]=median(subsq)
      
      if xhp+1 lt s0[0]-1 then im0[xhp+1,yhp]=median(im0[x1+1:x2+1,y1:y2])
      
    endfor

;
;wset, 0
;plot_image, bytscl(im0, 20, 200)

s=size(im0, /dim)

im1=im0

medcol=medcol2[*,im]

; first, add back on medcol2 values
for i=0, n_elements(medcol)-1 do im1[i,*]=im0[i,*]+medcol[i]

;wset, 1
;plot_image, bytscl(im1, 20, 200)

psf=modelpsf[*,*,im]

; now calculate new values, avoiding the target PSF!!
ipsf=where(psf ne 0, npsf)
i2d=array_indices(psf, ipsf)

; check location of target pixels
;oplot, i2d[0,*], i2d[1,*], color=cgcolor('purple'), psym=2

im1=float(im1)

; make PSF pixels NaN so median calculation treats these as missing data
for pix=0, npsf-1 do im1[i2d[0,pix], i2d[1,pix]]=!Values.F_NAN

; Update medcol3 with new values
for i=0, n_elements(medcol)-1 do medcol3[i,im]=median(im1[i,*])

; make a plot
fileout='~/BRITE/TESTSETS/werner4lc/plots/medcol_example_'+fname+'.ps'
chk=file_search(fileout, count=nchk)
if nchk eq 0 then medcol_plot, fileout, medcol, medcol3[*,im], im1, i2d, fname


median_diff[im]=total(medcol)-total(medcol3[*,im])

; now subtract new values from image
for i=0, n_elements(medcol)-1 do im1[i,*]=im0[i,*]+medcol2[i,im]
for i=0, n_elements(medcol)-1 do im1[i,*]=im1[i,*]-medcol3[i,im]

;;Check for negative pixels
;neg=where(im1 lt 0, nneg)
;
;if nneg gt 0 then begin
;
;  ineg=array_indices(im1, neg)
;
; ;replace negative pixels with zero
;  im1[ineg[0,*],ineg[1,*]]=0
;
;endif

im0=im1

; update data1 with new im0
data1[*,*,im]=im0

next_im:
endfor  ; end loop over images

; get new models from the updated images
for im=0, nimg-1 do begin
  
  if flag[im] ne 2 then continue
  
  im0=data1[*,*,im]
  
  ; check for negative pixels and make these 0
  negs=where(im0 lt 0, nneg)
  if nneg gt 0 then begin
    xneg=(array_indices(im0, negs))[0,*]
    yneg=(array_indices(im0, negs))[1,*]
    
    for nn=0, nneg-1 do im0[xneg[nn],yneg[nn]]=0
  endif
  
  negs=where(im0 lt 0, nneg)
  if nneg gt 0 then stop
  
  s=size(im0, /dim)
  
  ; add a border to the image
  im2=lonarr(s[0]+2,s[1]+2)
  im2[1,1]=im0
  
  ; use label region to determine number of illuminated regions and number of PSF pixels
  r2=label_region(im2 ge thr, /all_neighbors)
  r2=r2[1:s[0],1:s[1]]                         ; trim off border

  hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)

  ; sort hist2 in decreasing order
  sort2=reverse(sort(hist2))
  
  if n_elements(hist2) eq 1 then begin  ; no target OR bad
  flag[im]=0
  goto, next_im
  endif

; now get 2nd biggest group of pixels i.e. PSF pixels
  ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
  i2=array_indices(im0, ind2)

  xi=reform(i2[0,*])
  yi=reform(i2[1,*])

  ; update borders of PSF
  psf_loc[*,im]=[min(xi),max(xi),min(yi),max(yi)]

  npix_psf[im]=hist2[sort2[1]]

  ; record PSF for model
  for pix=0, npix_psf[im]-1 do newmodel[xi[pix],yi[pix],im]=data1[xi[pix],yi[pix],im]
  
;  ; check it
;  wset, 0
;  plot_image, bytscl(im0, 20, 200)
;  oplot, xi, yi, color=cgcolor('purple'), psym=2
;  stop
;  wset, 0
;  plot_image, bytscl(modelpsf[*,*,im], 20, 200)
;  oplot, xi, yi, color=cgcolor('purple'), psym=2
;  
;  wset, 1
;  plot_image, bytscl(newmodel[*,*,im], 20, 200)
;  oplot, xi, yi, color=cgcolor('purple'), psym=2
;  stop
  
endfor

; calculate average median column difference from median_diff
avg_med_diff=median(median_diff[where(median_diff ne 0)])
print, avg_med_diff

modelpsf=newmodel

; re-save .save file with xy_com added

save, filename=outdir+fname+'.sav', roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
  simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1, flag, $
  psf_loc, npix_psf, hp_loc, nhp, modelpsf, medcol3, avg_med_diff
  
endfor ; end loop over file

print, 'end of program'
end

