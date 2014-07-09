pro sig_mag_meet

; make scatter versus mag plots for meeting .ppt
; 
Compile_opt idl2

indir='/Users/gemmawhittaker/BRITE/data/UB/testing/nostack/p4_cust/'

outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/

filesin=file_search(indir+'*.sav', count=nsav)

; define variables
flx_scat=fltarr(nsav)
res=fltarr(nsav)
vmags=fltarr(nsav)

for i=0, nsav-1 do begin
  
  restore, filesin[i]
;  stop
  vmags[i]=vmag
  
  flux2=flux/robust_mean(flux,2)
  
  flx_scat[i]=robust_sigma(flux2)
  
  res[i]=robust_mean(mean_res,2)
   
endfor

fileout0=outdir+'gem_scat_mag.ps'
plotsym, 0, /fill, 0.5
ps_on, fileout0, xsize=15, ysize=15
plot, vmags, flx_scat, xtitle='V-Magnitude', ytitle='Scatter', charsize=0.7, psym=8
ps_off

;fileout1=outdir+'gem_res_scat.ps'
;ps_on, fileout1, xsize=15, ysize=15
;plot, res, flx_scat, xtitle='Residual', ytitle='Scatter', charsize=0.7, psym=2
;ps_off

stop
end 