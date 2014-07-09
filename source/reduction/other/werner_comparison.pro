pro werner_comparison

; produce comparison results for old data versus new data
; 
Compile_opt idl2

!p.background=cgcolor('white')
plotsym, 0, /fill, 0.7

cols=[cgcolor('black'), cgcolor('purple')]

indir1='/Users/gemmawhittaker/BRITE/data/UB/p4_1/'
indir2='/Users/gemmawhittaker/BRITE/data/UB/p4_3/'

outdir='/Users/gemmawhittaker/BRITE/reports/myreports/'
 
filesin1=file_search(indir1+'*.sav', count=nf1)
filesin2=file_search(indir2+'*.sav', count=nf2)

mf1=fltarr(nf1)
mf2=fltarr(nf2)

mags=fltarr(nf1)

scat1=fltarr(nf1)
scat2=fltarr(nf2)

for i=0, nf1-1 do begin
  
  restore, filesin1[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
  f1=flux
  t1=jd
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
    
    fmag=-2.5*alog10(f_temp)
    
    f_norm=fmag/robust_mean(fmag,2)
    
    fscat1[j]=stddev(f_norm)
    
    mean_mag[j]=robust_mean(fmag,2)
    
    stop
    
  endfor
  
  ; calulculate average scatter across all orbits and mean magnitude
  mf1[i]=robust_mean(mean_mag,2)
  scat1[i]=robust_mean(fscat1,2)
  
  
  restore, filesin2[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
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
    
    fmag=-2.5*alog10(f_temp)
    
    f_norm=fmag/robust_mean(fmag,2)
    
    fscat2[j]=stddev(f_norm)
    
    mean_mag[j]=robust_mean(fmag,2)
    
  endfor
  
  ; calulculate average scatter across all orbits and mean magnitude
  mf2[i]=robust_mean(mean_mag,2)
  scat2[i]=robust_mean(fscat2,2)
  
endfor

; make plots
fileout=outdir+'v1_v2_stddev.ps'
!p.multi=[0,1,2,0,0]
ps_on, fileout, xsize=15, ysize=25
plot, mags, mf1, psym=8, color=cgcolor('black'), xtitle='V-mag', ytitle='Instrumental Mag', charsize=0.7, $
  title='UB Photometry - 1s exp - V1 vs V2'
oplot, mags, mf2, psym=8, color=cgcolor('purple')
al_legend, ['v1', 'v2'], psym=[8,8], colors=cols, charsize=0.7

plot, mags, scat1, psym=8, color=cgcolor('black'), xtitle='V-mag', ytitle='STDDEV', charsize=0.7
oplot, mags, scat2, psym=8, color=cgcolor('purple')
al_legend, ['v1', 'v2'], psym=[8,8], colors=cols, charsize=0.7

ps_off

spawn, 'open '+fileout+' &'


stop
print, 'End of Program'
end