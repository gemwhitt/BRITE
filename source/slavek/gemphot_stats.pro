pro gemphot_stats

Compile_opt idl2

; program to read in results from p3_p4_gem3 and plot stats:
; 
; flux - x4 - scatter 
; number of pixels in PSF - x4 - scatter
; max DN in brightest pixel - scatter
; max SNR - scatter
; Vmag
; Instrumental mag

indir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot4/'
outdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot4/stats/'

filesin=file_search(indir+'HD*.sav', count=nfiles)

opt=long(0)

fileout=outdir+'1s_stats_'+strtrim(opt,2)+'.txt'

for i=0, nfiles-1 do begin
  
  restore, filesin[i] ;vmag, nsubpix, max_snr, max_dn, xy_psf, flux, jd, nsat, $
                      ;roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
                      ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend
                      
  hdname=file_basename(filesin[i], '_p4_gem1.sav')
  
  nfrm=n_elements(jd)
  
  ; do scatter etc per orbit - then get average
  jd2=jd[1:nfrm-1]
  jdiff=jd2-jd
  gap=where(jdiff gt 0.015, ngap)
  gap=[0,gap, nfrm-1]
  
  imag=-2.5*alog10(flux)
  
  ; variables.... ;;;;;;;;;;;;;;;;;
  mean_npix=fltarr(ngap+1)
  sig_npix=fltarr(ngap+1)
  
  mean_snr=fltarr(ngap+1)
  sig_snr=fltarr(ngap+1)
  
  mean_dn=fltarr(ngap+1)
  sig_dn=fltarr(ngap+1)
  
  mean_flux=fltarr(ngap+1)
  sig_flux=fltarr(ngap+1)
  
  ;vmag
  ;nsat
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  for j=0, ngap do begin
    
    if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+(gap[j]+1)
    
    mean_npix[j]=robust_mean(nsubpix[opt,iloc],2)
    sig_npix[j]=robust_sigma(nsubpix[opt,iloc]/mean_npix[j])
    
    mean_snr[j]=robust_mean(max_snr[iloc],2)
    sig_snr[j]=robust_sigma(max_snr[iloc]/mean_snr[j])
    
    mean_dn[j]=robust_mean(max_dn[iloc],2)
    sig_dn[j]=robust_sigma(max_dn[iloc]/mean_dn[j])
    
    mean_flux[j]=robust_mean(flux[opt,iloc],2)
    sig_flux[j]=robust_sigma(flux[opt,iloc]/mean_flux[j])
 
;if j ge 15 then stop    
  endfor
 
  ; get robust means of all stats
  mean_npix=robust_mean(mean_npix,2)
  sig_npix=robust_mean(sig_npix,2)
  
  mean_snr=robust_mean(mean_snr,2)
  sig_snr=robust_mean(sig_snr,2)
  
  mean_dn=robust_mean(mean_dn,2)
  sig_dn=robust_mean(sig_dn,2)
  
  mean_flux=robust_mean(mean_flux,2)
  sig_flux=robust_mean(sig_flux,2)
  
  nsat=robust_mean(nsat,2)
  
  vmag=vmag
;  stop
  ; save variables
  openw, lun, fileout, /get_lun, /append
  printf, lun, hdname, vmag, nsat, mean_npix, sig_npix, mean_snr, sig_snr, mean_dn, sig_dn, mean_flux, sig_flux, $
    format='(a7, x, f7.3, x, i3, x, f, x, f, x, f, x, f, x, f, x, f, x, f, x, f)'
  free_lun, lun
  
 
endfor
;spawn, 'open '+fileout+' &'

stop

end
