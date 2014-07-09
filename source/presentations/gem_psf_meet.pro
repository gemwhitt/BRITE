pro gem_psf_meet

  ; do psf-fitting for p3_test data - analyse the results
  ;
  Compile_opt idl2
  
  timefile='~/Desktop/timefile_custom.dat'
  ;
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  ;
  toplot='y'
  ;
  ;
  indir='~/BRITE/data/UB/testing/nostack/p3/'
  filesin=file_search(indir+'*.sav', count=nfiles)
  
  ;outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/'
  outdir='/Users/gemmawhittaker/BRITE/data/UB/testing/nostack/p4_cust/'
  
  ; control parameters:
  rbin=8
  norb=1
  
  for gg=2, nfiles-1 do begin
  
    print, file_basename(filesin[gg],'.sav')+' - file '+strtrim(gg,2)
    
    restore, filesin[gg] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, sm_flx, avg_flx, expn, $
   ;   simbad_radec, simbad_mags, parlax, otype, sptype
   
   vpos=strpos(simbad_mags, 'V=')
   vmag=float(strmid(simbad_mags, vpos+2))
   
    
    dark_current=0. ; change this when dark current is known
    
    fname=file_basename(filesin[gg],'.sav')
    dashpos=strpos(fname, '_')
    hdname=strmid(fname, 0, dashpos)
    
    dat=data1 ; no effect to data1
    
    ; get total number of frames in this file
    nfrm=(size(dat, /dim))[2]
    
    flux=fltarr(nfrm)
    backgd=fltarr(nfrm)
    backer=fltarr(nfrm)
    mean_res=fltarr(nfrm)
    sig_res=fltarr(nfrm)
    
    ; get number of orbits
    jd2=jd[1:nfrm-1]
    timediff=jd2-jd
    
    starttime=systime()
    for img=0, nfrm-1 do begin  ;ngap-1 do begin
    
      it_pks=0
    
      print, img
      redo=0
      
      dat1=dat[*,*,img]
      
      ; rebin the image - 4x4, 8x8, 16x16...
      dat2=rebin_data(dat1,rbin)
      
      ndat2=n_elements(dat2[0,0,*])  ; number of binned pixels
      
      ; find approx center of PSF
      mthr=200
      spr=2*rbin
      
      get_peaks:
      pks=brite_findpeaks(dat2, minthresh=mthr, spread=spr)
      
      pks1=round(pks)
      
      ; determine the background - away from the PSF center
      get_bkgd, dat2,pks,rbin, bkgd,bk_err
      
      ; subtract the background measurement from the image
      dat2=dat2-bkgd
      
      backgd[img]=bkgd
      backer[img]=bk_err
      
      npks=n_elements(pks)/3.  ; number of bright peaks found in image - should be only 1.
      
      if npks gt 1 then begin
        ; if more than one peak then process accordingly....
        brightest=where(pks[2,*] eq max(pks[2,*]))
        pks=pks[*,brightest]
        pks1=pks
      endif  
      
      if  pks[0] eq -999.9 then begin
        it_pks=it_pks+1
        if it_pks ge 2 then stop
        spr=3*rbin
        goto, get_peaks
      endif
      
      ; make cutout - centered on PSF center
      if pks1[0]-(rbin*6) le 0 then x0=0 else x0=pks1[0]-(rbin*6)
      if pks1[0]+(rbin*6) ge rbin*32. then x1=(rbin*32)-1 else x1=pks1[0]+(rbin*6)
      if pks1[1]-(rbin*6) le 0 then y0=0 else y0=pks1[1]-(rbin*6)
      if pks1[1]+(rbin*6) ge rbin*32. then y1=(rbin*32)-1 else y1=pks1[1]+(rbin*6) 
      
      cutout=dat2[x0:x1,y0:y1]
      
 ;     plot_image, bytscl(cutout, 50, 5000)
      
;      stop
      
      redo_start:
      
      ncut=n_elements(cutout)
      
      ;reorder pixels by data number
      order1=reverse(sort(cutout))
      temp_arr=order1
      signal1=cutout[order1]
      
      count=-1
      snr1=0           ; temp removed
      ;snr1=fltarr(ncut)  ; extra added
      signal=fltarr(ncut)
      
      repeat begin
      
        count=count+1
        
        snr=snr1
        
        sig=total(signal1[0:count])
        
        signal[count]=sig
        
        ;snr1=sig/sqrt(count+1) ;*((bkgd + dark_current)^2))       ; temp removed
        snr1=sig/sqrt(sig + (count+1)*(50)^2)       ; test
        ;snr1[count]=sig/sqrt(sig + (count+1)*((bkgd + dark_current)^2))   ; extra edded
        
        temp_arr[count]=-9999
        
      endrep until snr1 lt snr      ; temp removed
      ;endrep until count eq ncut-2   ; extra added
      
      if count ge ncut-100 then begin
        stop
 ;       redo=redo+1
 ;       cutout=(dat2[pks1[0]-(rbin*(6+redo)):pks1[0]+(rbin*(6+redo)),pks1[1]-(rbin*(6+redo)):pks1[1]+(rbin*(6+redo))])
 ;       goto, redo_start
      endif
      
      
      ; get mean and scatter of residuals in the ROI....
      temp_dat2=dat2
      subsq=[x0,x1,y0,y1]
      temp_dat2[subsq[0:1],subsq[2:3]]=-9999
      
      bb=where(temp_dat2 ne -9999)
      ;print, robust_mean(dat2[bb],2)
      
      respix=[signal1[count+1:ncut-1],dat2[where(temp_dat2 ne -9999)]]
      
      loc_2d=array_indices(cutout, order1)
      
      xx=where(temp_arr eq -9999, complement=yy)
      
     ; plot_image, bytscl(cutout, 50, 500)
     ; oplot, loc_2d[0,xx], loc_2d[1,xx], psym=2, color=cgcolor('orange')
     ; oplot, loc_2d[0,yy], loc_2d[1,yy], psym=2, color=cgcolor('purple')
      
           ;  stop
      
      mean_res[img]=robust_mean(respix,2)
      sig_res[img]=robust_sigma(respix)
      
      flux[img]=sig
      
      maxloc=where(snr1 eq max(snr1), nmax)
      if nmax gt 1 then maxloc=maxloc[nmax-1]
      ;       stop
      
      ;fileout1=outdir+'gem_psf.ps'
      ;  ps_on, fileout1, xsize=15, ysize=11
      ;  plot_image, bytscl(cutout, 50, 500), title=hdname, charsize=0.7
        ;oplot, loc_2d[0,xx], loc_2d[1,xx], psym=4, color=cgcolor('orange'), symsize=2
      ;  oplot, loc_2d[0,0:maxloc], loc_2d[1,0:maxloc], psym=4, color=cgcolor('orange'), symsize=2
        ;plot, snr1, color=cgcolor('black'), xtitle='Number sub-pixels contained in PSF', ytitle='SNR', charsize=0.7
      ;  ps_off
        
;fileout2=outdir+'gem_snr.ps'
;!x.thick=2
;!y.thick=2
       ;ps_on, fileout2, xsize=15, ysize=15
       ; plot, snr1, color=cgcolor('black'), xtitle='Number sub-pixels contained in PSF', ytitle='SNR', charsize=0.7, $
       ;   thick=3
       ; ps_off
        
        
;        stop
      
      
      ;print, 'max at ',where(snr1 eq max(snr1))
      ;print, count
      
      ;stop
      
    endfor ; END LOOP OVER THIS SUB-SECTION OF ORBITS - IS USING
    endtime=systime()
    
   ; stop
    
    ;openw, lun, timefile, /append, /get_lun
    ;printf, lun, hdname, rbin, starttime, endtime, format='(a,2x,i,2x,a,2x,a)'
    ;free_lun, lun
    
    bkgd=backgd
    bk_err=backer
    
     
    save, filename=outdir+hdname+'_p4_test.sav', flux, jd, bkgd, bk_err, mean_res, sig_res, vmag
    
endfor  ; end loop over file
  stop
end