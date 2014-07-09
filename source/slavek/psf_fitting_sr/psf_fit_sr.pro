pro psf_fit_sr

; do psf-fitting for p3 data - analyse the results
; 
Compile_opt idl2
;
; FOR PLOTTING;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')
;
todisplay='n'  ; see Step 5
;
;
indir='~/BRITE/data/UB/p3/'
filesin=file_search(indir+'Orion-CF1-2*.sav', count=nfiles)

outdir='~/BRITE/data/UB/p4_sr1/'
norb=1  ; number of orbits (45 frames per orbit) to use to get average PSF

runtime=fltarr(nfiles)

for gg=0, nfiles-1 do begin
  
  stime=systime(1)
  
  print, file_basename(filesin[gg],'.sav')+' - file '+strtrim(gg,2)

  restore, filesin[gg] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
  jd=jd-jd[0]
  keep=where(jd gt 0.)
  
  roi_name=roi_name[keep]
  exp_num=exp_num[keep]
  jd=jd[keep]
  ccd_temp=ccd_temp[keep]
  p2_trend=p2_trend[keep]
  data1=data1[*,*,keep]
  
  fname=file_basename(filesin[gg],'.sav')
  dashpos=strpos(fname, '_')
  hdname=strmid(fname, 0, dashpos)
  
  ; get Vmag
  vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
 
  dat=data1 ; no effect to data1
  
  ; get total number of frames in this file
  nfrm=(size(dat, /dim))[2]

  temp_arr=indgen(nfrm)
  
    
  ; get number of orbits
  jd2=jd[1:nfrm-1]
  timediff=jd2-jd
  
  new_orbit=where(timediff gt 0.015, ngap)  ; there will be ngap+1 loops
  new_orbit=[new_orbit,nfrm-1]
  temp_orbit=new_orbit
  for rr=0+(norb-1), ngap, norb do temp_orbit[rr]=-999
  keep=where(temp_orbit eq -999)
  new_orbit=new_orbit[keep]
  ngap=n_elements(new_orbit)-1
  
  rbin=8.
  
  xy_psf=fltarr(2,nfrm)
  
  for orbit=-1, ngap-1 do begin
      
    starttime=systime()
        
    if orbit eq -1 then iloc=indgen(new_orbit[0]+1) else iloc=indgen(new_orbit[orbit+1]-new_orbit[orbit])+new_orbit[orbit]+1
    
;    if orbit eq 38 then stop

    dat1=dat[*,*,iloc]
    
    day=jd[iloc]
   
    
    nfrm=(size(dat1, /dim))[2]
 
    ; find approx center of PSF
    mthr=400.
    spr=2
    x = fltarr(32,32)                   ; average image
    m = nfrm                                      ; how many images
    for i=0,m-1 do x = x + dat1[*,*,i]            ; adding images
    x = x/m                                       ; average image
    pks=brite_findpeaks(x, minthresh=mthr, spread=spr)
    
    npks=n_elements(pks)/3.
    
    if pks[0] eq -999.9 then temp_arr[iloc]=-9999
    if pks[0] eq -999.9 then continue
   
    if npks gt 1 then begin
      brightest=where(pks[2,*] eq max(pks[2,*]))
      pks=pks[*,brightest]
    endif
    
    pks2=pks
    pks2[0]=pks[0]*rbin
    pks2[1]=pks[1]*rbin
 
    ; rebin the images - 4x4, 8x8, 16x16...
    dat2=rebin_data(dat1,rbin)             ; ORIGINAL!
          
    ; Step 3: Use a Gaussian as the 0-th approximation for the PSF...
    cen=8000.
    sig=2.
    psf0 = first_approx(cen,sig,rbin,pks2)
       
    ; Step 4: Determine the shifts in small pixels from the image centre for individual images assuming a given PSF.
    ; note that these shifts are given in "small pixels", so shx/8. is the shift in big pixels
    get_shifts, dat2,psf0,rbin, shx1,shy1

    ; Step 5: Optional step - useful to see what is happening, e.g. one can check how well the centres are found.
    if todisplay eq 'y' then $
      for i=0, 10 do begin
        shx=shx1[i]
        shy=shy1[i]
        
        
        window, 0, xpos=1500, ypos=200
        plot_image, x
        oplot, [pks[0]], [pks[1]], psym=2, color=cgcolor('purple')
        
        window, 1, xpos=2300, ypos=200
        plot_image, dat2[*,*,i]
        oplot, [pks2[0]+shx[0]], [pks2[1]+shy[0]], psym=2, color=cgcolor('purple')
        
        window, 2, xpos=2000, ypos=-400
        plot_image, psf0      
        oplot, [pks2[0]], [pks2[1]], psym=2, color=cgcolor('purple')
stop
    endfor
    
    count=-1
    repeat begin
  
      count=count+1
      shx=shx1
      shy=shy1
    
      ; Step 6: Using the shx, shy, adjust individual images to the centre. The determine the improved 2-D PSF.
      ; Note the background is subtracted (for use as an PSF), but is given by "bck".
      psf1 = get_psf(dat2,shx,shy,rbin, dud,bck,w)  ; w is 2d-indices of background pixels
      
      xw=(array_indices(findgen(rbin*32,rbin*32), w))[0,*]
      yw=(array_indices(findgen(rbin*32,rbin*32), w))[1,*]
      
      ;window, 2, xpos=2000, ypos=-300
      ;for i=0, n_elements(shx)-1 do begin
      ;plot_image, dat2[*,*,i]
      ;oplot, [shx[i]]+pks[0], [shy[i]]+pks[1], color=cgcolor('purple'), psym=2
      ;oplot, xw, yw, color=cgcolor('orange'), psym=1
      ;wait, 1
      ;endfor

      get_shifts, dat2,psf1,rbin, shx1,shy1
 ;stop     
      ; compute the residuals between shifts of old fit and shifts of new fit
      shx_res=abs(shx1-shx)
      shy_res=abs(shy1-shy)
    
      ; determine if there are any shifts and if so, repeat
      rs1=where(shx_res gt 0, nrs1)
      rs2=where(shy_res gt 0, nrs2)
    
      if nrs1 gt 0 OR nrs2 gt 0 then cycle=1 else cycle=0
    
    endrep until cycle eq 0  ; step 7 is to iterate until convergence
    
    ; count = number of iterations
    ;print, 'Number of iterations was '+strtrim(count,2)
  
   ; optional display....to see PSF
   if todisplay eq 'y' then begin
    psf_disp,psf1,'PSF-1',0.,0.,rbin
   ; stop
   endif
;  stop
    ; Step 8: The flux & background determination by linear fits of the experimental PSF.
    psf3=psf1;[pks[0]-(5*rbin):pks[0]+(5*rbin),pks[1]-(5*rbin):pks[1]+(5*rbin)]
    shx2=shx
    shy2=shy
    dud2=dud
    
    ; save PSF center locations
    
;    for i=0, 10 do begin
;     
;      window, 1, xpos=2300, ypos=200
;      plot_image, dat2[*,*,i]
;      oplot, [pks2[0]+shx[i]], [pks2[1]+shy[i]], psym=2, color=cgcolor('purple')
      
;      stop
;    endfor
   
    xy_psf[0,iloc]=(pks2[0]+shx)/8.
    xy_psf[1,iloc]=(pks2[1]+shy)/8.
   
    ; plot_image, psf3
   ;stop
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; this will be end of loop for finding the psf
    ;
    wgt=sqrt(dat1)
 
    fit = B2_im8(psf3,dat2,wgt,shx2,shy2,dud2)  ; result of function: coefficients of fit as 5-column array
;          [flux, background, error flux, err-backgr, err-per-pixel]
      
    ;fastphot, reform(dat2[*,*,i]), psf3, noise_map, xloc, yloc, flux, map_clean=map_clean, eflux=eflux, mu=mu
 
    if orbit eq -1 then begin
      flux=reform(fit[0,*])
      ferr=reform(fit[2,*])
      bkgd=reform(fit[1,*])
      berr=reform(fit[3,*])
      pixerr=reform(fit[4,*])
    endif else begin
      flux=[flux,reform(fit[0,*])]
      ferr=[ferr,reform(fit[2,*])]
      bkgd=[bkgd,reform(fit[1,*])]
      berr=[berr,reform(fit[3,*])]
      pixerr=[pixerr,reform(fit[4,*])]
    endelse
   
; stop
   endfor ; END LOOP OVER THIS SUN-SECTION OF ORBITS - IS USING
   
   rej=where(temp_arr eq -9999, nrej, complement=keep)
   jd=jd[keep]
   ccd_temp=ccd_temp[keep]
   p2_trend=p2_trend[keep]
   xy_psf=xy_psf[*,keep]
   data1=data1[*,*,keep]
   
   
   ; save results
   fileout=outdir+hdname+'_p4_sr.sav'
   save, filename=fileout, vmag, flux, ferr, bkgd, berr, pixerr, jd, $
     roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
     simbad_mags, parlax, otype, sptype, p2_trend, med_trend, xy_psf, data1
     
     etime=systime(1)
    
;stop   
runtime=etime-stime
endfor  ; end loop over file

print, 'total run time is '+strtrim(total(runtime)/60.,2)+' mins'


;stop
end