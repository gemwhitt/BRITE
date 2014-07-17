pro cf_scatter

; progam to compare light curves from Gemma and Slavek photometry....
; 
Compile_opt idl2

!p.background=cgcolor('white')
plotsym, 0, /fill, 0.8

ingem='~/BRITE/data/UB/p4_gemphot1/'
inslr='~/BRITE/data/UB/p4_sr1/'

fgem=file_search(ingem+'*.sav', count=nfiles)
fslv=file_search(inslr+'*.sav', count=nfiles)

cols1=[cgcolor('purple'), cgcolor('sea green')]
cols2=['purple', 'sea green']
  
outdir='/Users/gemmawhittaker/BRITE/reports/myreports/'
  
mags=fltarr(nfiles)
  
mf1=fltarr(nfiles)
mf2=fltarr(nfiles)
scat1=fltarr(nfiles)
scat2=fltarr(nfiles)
  
for i=0, nfiles-1 do begin
  
restore, fgem[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
t1=jd

fscat=fltarr(4)
for j=0, 3 do begin
  meanf=robust_mean(flux[j,*],2)
  normf=flux[j,*]/meanf
  fscat[j]=robust_sigma(normf)
endfor

xx=(where(fscat eq min(fscat)))[0]
f1=flux[xx,*]

mags[i]=vmag
n1=n_elements(t1)
    
  t11=t1[1:n1-1]
  dif1=t11-t1
  gap=where(dif1 gt 0.015, ngap1)
  gap=[0,gap,n1-1]
    
  fscat1=fltarr(ngap1+1)
  mean_mag=fltarr(ngap1+1)
    
    for j=0, ngap1 do begin
    
      if j eq 0 then f_temp=f1[0:gap[j+1]] else f_temp=f1[gap[j]+1:gap[j+1]]
      
      fmag=f_temp
      
      f_norm=fmag/robust_mean(fmag,2)
      
      fscat1[j]=stddev(f_norm)
      
      mean_mag[j]=robust_mean(fmag,2)
      
    endfor
    
    ; calulculate average scatter across all orbits and mean magnitude
    mf1[i]=robust_mean(mean_mag,2)
    scat1[i]=robust_mean(fscat1,2)
    
    
    restore, fslv[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
    f2=flux
    t2=jd
    n2=n_elements(t2)
    
    t22=t2[1:n2-1]
    dif2=t22-t2
    gap=where(dif2 gt 0.015, ngap2)
    gap=[0,gap,n2-1]
    
    fscat2=fltarr(ngap2+1)
    mean_mag=fltarr(ngap2+1)
    
    for j=0, ngap2 do begin
    
      if j eq 0 then f_temp=f2[0:gap[j+1]] else f_temp=f2[gap[j]+1:gap[j+1]]
      
      fmag=f_temp
      
      f_norm=fmag/robust_mean(fmag,2)
      
      fscat2[j]=stddev(f_norm)
      
      mean_mag[j]=robust_mean(fmag,2)
      
    endfor
    
    ; calulculate average scatter across all orbits and mean magnitude
    mf2[i]=robust_mean(mean_mag,2)
    scat2[i]=robust_mean(fscat2,2)
    
  endfor
  
  ; make plots
  fileout=outdir+'gem_slv_stddev.ps'
 ; !p.multi=[0,1,2,0,0]
  ps_on, fileout, xsize=15, ysize=15
;  plot, mags, mf1, psym=8, color=cgcolor('black'), xtitle='V-mag', ytitle='Instrumental Mag', charsize=0.7, $
;    title='UB Photometry - 1s exp - V1 vs V2'
;  oplot, mags, mf2, psym=8, color=cgcolor('purple')
;  al_legend, ['v1', 'v2'], psym=[8,8], colors=cols, charsize=0.7
  
  plot, mags, scat2, psym=8, color=cgcolor('black'), xtitle='V-mag', ytitle='STDDEV', charsize=0.7, /nodata
  oplot, mags, scat1, psym=8, color=cols1[0]
  oplot, mags, scat2, psym=8, color=cols1[1]
  al_legend, ['gem', 'sr'], psym=[8,8], colors=cols2, charsize=0.7
  
  ps_off
  
  spawn, 'open '+fileout+' &'
  
  
  stop
  print, 'End of Program'
end