pro save_psf_info

; Main program to run get_psf_info.pro and save results into a .txt or .dat file 
; to be called by p4_gem7.pro
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.8, /fill
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

nstk=5
!p.multi=0

indir='~/BRITE/data/UB/p2/CENTAURUS/'+strtrim(nstk,2)+'stk/'

filesin=file_search(indir+'*.sav', count=nfiles)

fname=file_basename(filesin, '_p2_'+strtrim(nstk,2)+'.sav')

outdir='/Users/gemmawhittaker/BRITE/data/UB/p4/CENTAURUS/'+strtrim(nstk,2)+'stk/psf/'

for b=0, nfiles-1 do begin

  print, fname[b]
  
  restore, filesin[b] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
  ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
  
  bminv=bmag-vmag
  
  hdname=roi_name[0]
  
  jd1=jd-jd[0]
  keep=where(jd1 ge 0.0)
  
  bad=where(jd1 lt 0.0, nbad)
  if nbad gt 0 then stop
  
  dat=data1 ; no effect to data1 OR times 10. for 0.1s exp
  
  ;get total number of frames in this file
  nfrm=(size(dat, /dim))[2]
  
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  ; calculate cadence of observations
  cadence=robust_mean(jdiff,2)  ; in days
  cadence_m=cadence*24.*60.       ; in minutes
  cadence_s=cadence*24.*3600.       ; in seconds
  n_img_per_orbit=fix(15./cadence)
  
  max_dn=lonarr(nfrm)
  
  for i=0, nfrm-1 do max_dn[i]=max(dat[*,*,i])
  
  xy_psf=fltarr(2,nfrm)
  npix=fltarr(3)  ; value plux sigma
  e_npix=fltarr(3) ; value plux sigma
  fl=dblarr(3)  ; value plux sigma
  e_fl=dblarr(3) ; value plux sigma
  fwhm=fltarr(2)
    
  ; call external program - get PSF centers, average flux per orbit and number of pixels in PSF
  get_psf_info, dat,gap,ngap,nfrm,outdir,hdname, xy_psf,npix,e_npix,fl,e_fl,fwhm
  
  ;get total number of frames in this file
  nfrm=(size(dat, /dim))[2]
  
  ; save results
  fileout=outdir+fname[b]+'_psf.txt'
  
  openw, lun, fileout, /get_lun
  printf, lun, 'FWHM',fwhm[0],fwhm[1], format='(a15,x,f7.3,x,f7.3)'
  printf, lun, 'NUM_PIX_IN_AP',npix[0],npix[1],npix[2], format='(a15,x,i7,x,i7,x,i7)'
  printf, lun, 'E_NUM_PIX',e_npix[0],e_npix[1],e_npix[2], format='(a15,x,f7.4,x,f7.4,x,f7.4)'
  printf, lun, 'AVG_FLUX',fl[0],fl[1],fl[2], format='(a15,x,d14.1,x,d14.1,x,d14.1)'
  printf, lun, 'E_FLUX',e_fl[0],e_fl[1],e_fl[2], format='(a15,x,f7.4,x,f7.4,x,f7.4)'
  printf, lun, 'VMAG',vmag, format='(a10,x,f7.3)'
  printf, lun, 'B-V_MAG',bminv, format='(a15,x,f7.3)'
  printf, lun, ''
  for i=0, nfrm-1 do printf, lun, max_dn[i], xy_psf[0,i], xy_psf[1,i], format='(i10, x, f7.2, x, f7.2)'
  free_lun, lun
  
  
endfor

print, 'end of program'
end
  
