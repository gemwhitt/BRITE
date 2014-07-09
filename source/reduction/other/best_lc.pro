pro best_lc

  Compile_opt idl2
  
  ; choose best LC based on gemphot stats
  
  indir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot3/stats/'
  outdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot3/stats/'
  
  filesin=file_search(indir+'1s*.txt', count=nfiles)
  
  opt='stk' ;'stk'
  
  sig_lc=fltarr(15,4)
  npix=fltarr(15,4)
  sig_pix=fltarr(15,4)
  lc_adu=fltarr(15,4)
    
  for i=0, nfiles-1 do begin
  
    if opt eq 'sin' then $
      readcol, filesin[i], hdname, vmag, nsat, mean_npix, sig_npix, mean_snr, sig_snr, mean_dn, sig_dn, mean_flux, sig_flux, $
      format='a, f, i, f, f, f, f, f, f, f, f', numline=15 else $
      readcol, filesin[i], hdname, vmag, nsat, mean_npix, sig_npix, mean_snr, sig_snr, mean_dn, sig_dn, mean_flux, sig_flux, $
      format='a, f, i, f, f, f, f, f, f, f, f', skipline=15
   
   
   sig_lc[*,i]=sig_flux 
   npix[*,i]=mean_npix
   sig_pix[*,i]=sig_npix
   lc_adu[*,i]=mean_flux
        
  endfor

  best=intarr(15)
  
  for i=0, 14 do best[i]=where(sig_lc[i,*] eq (min(sig_lc[i,*]))[0])
  
  fileout=outdir+'best_lcs.txt'
  openw, lun, fileout, /get_lun, /append
  for i=0, 14 do printf, lun, hdname[i], vmag[i], best[i], nsat[i], lc_adu[i,best[i]], sig_lc[i,best[i]], $
    npix[i, best[i]], sig_pix[i, best[i]], mean_dn[i], sig_dn[i], $
    format='(a7, x, f, x, i1, x, i3, x, f, x, f, x, f, x, f, x, f, x, f)'
    free_lun, lun
    
    spawn, 'open '+fileout+' &'
  
end
