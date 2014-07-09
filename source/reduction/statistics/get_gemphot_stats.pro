pro get_gemphot_stats

; program to measure statistics from results of p4_gemphotXX - then plot with plot_gemphpt_stats.pro
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.8, /fill
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
  
indir='~/BRITE/data/UB/p4_gemphot3b/'
filesin=file_search(indir+'Orion-CF1-3*.sav', count=nfiles)
;filesin=file_search(indir+'HD*.sav', count=nfiles)
fname=file_basename(filesin, '_p4.sav')
  
outdir=indir+'stats/'

opt=long(3)

fileout=outdir+'1s_stats_'+strtrim(opt,2)+'.txt'

meanflux=fltarr(15)
eflux=fltarr(15)
pix=intarr(15)
epix=fltarr(15)
maxdn=fltarr(15)
emaxdn=fltarr(15)
vmags=fltarr(15)
name=strarr(15)
nsat1=intarr(15)
    
for bb=0, nfiles-1 do begin
  
    restore, filesin[bb] ;vmag, nsubpix, max_snr, max_dn, xy_psf, flux, jd, nsat, $
                         ;roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
                         ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend, exp_ttl, exp_time
                         
    ;wset, 1
    ;plot, jd, flux[0,*], color=cgcolor('black'), psym=2, /ynozero
                         
                        
    name[bb]=roi_name[0]
    vmags[bb]=vmag
    
    nsat1[bb]=median(nsat)
    
    nfrm=n_elements(jd)
    
    jd1=jd
    jd2=jd1[1:nfrm-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap)
    gap=[0,gap,nfrm-1]
    
    sub_flux=fltarr(ngap+1)
    sub_eflux=fltarr(ngap+1)
    sub_pix=intarr(ngap+1)
    sub_epix=fltarr(ngap+1)
    sub_maxdn=fltarr(ngap+1)
    sub_emaxdn=fltarr(ngap+1)
    
    for i=0, ngap do begin
    
      if i eq 0 then iloc=indgen(gap[i+1]+1) else iloc=indgen(gap[i+1]-gap[i])+gap[i]+1
          
      sub_flux[i]=median(flux[opt,iloc])
      sub_eflux[i]=robust_sigma(flux[opt,iloc]/sub_flux[i])
      sub_pix[i]=median(nsubpix[opt,iloc])
      sub_epix[i]=robust_sigma(nsubpix[opt,iloc]/sub_pix[i])
      sub_maxdn[i]=median(max_dn[iloc])
      sub_emaxdn[i]=robust_sigma(max_dn[iloc]/sub_maxdn[i])
          
    endfor  ; end loop over gap
    
    meanflux[bb]=median(sub_flux)
    eflux[bb]=median(sub_eflux)
    pix[bb]=median(sub_pix)
    epix[bb]=median(sub_epix)
    maxdn[bb]=median(sub_maxdn)
    emaxdn[bb]=median(sub_emaxdn)
    
  endfor  ; end loop over file
  
  ;write out data
  openw, lun, fileout, /get_lun
  for i=0, 14 do printf, lun, name[i], vmags[i], nsat1[i], pix[i], epix[i], maxdn[i], emaxdn[i], meanflux[i], eflux[i], $
    format='(a10, x, f, x, i, x, i, x, f, x, f, x, f, x, f, x, f)'
    free_lun, lun
    
;    hdname, vmag, nsat, mean_npix, sig_npix, mean_snr, sig_snr, mean_dn, sig_dn, mean_flux, sig_flux, $
;      format='a,f,i,f,f,f,f,f,f,f,f'

  print, 'end of program'
;  stop
end


