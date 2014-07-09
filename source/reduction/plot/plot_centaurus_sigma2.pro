pro plot_centaurus_sigma2

; 2nd effort
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.8, /fill

nstk=['1','5','10']

indir='~/BRITE/data/UB/p4/CENTAURUS/'+nstk+'stk/lc_txt/'
ndir=n_elements(indir)

outdir='~/Desktop/'

fileout=outdir+'lc_sigmas.txt'
plotout=outdir+'lc_sigma.pdf'

cols1=cgcolor(['blue', 'green', 'red'])

mean_mag=fltarr(ndir,30)
dmag=fltarr(ndir,30)
vmags=fltarr(ndir,30)

for d=0, ndir-1 do begin

  filesin=file_search(indir[d]+'*.txt', count=nf)
  
  for i=0, nf-1 do begin
  
    readcol, filesin[i], skipline=14, jd, mag, format='d,f,x,x,x,x,x,x,x'
    
    readcol, filesin[i], skipline=6, numline=1, vmag, format='a'
    vmag=strtrim(vmag,2)
    eqpos=strpos(vmag, '=')
    vmag=float(strmid(vmag, eqpos+1))
    
    vmags[d,i]=vmag
    
    n=n_elements(jd)
    
    jd1=jd-jd[0]
    jd2=jd1[1:n-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap)
    gap=[0,gap, n-1]
    
    sigmag=fltarr(ngap+1)
    
    for j=0, ngap do begin
      
      if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+gap[j]+1
      
      fl=mag[iloc]
      medfl=robust_mean(fl,2)
      
      sigmag[j]=robust_sigma(fl/medfl)
      
    endfor
    
    mean_mag[d,i]=robust_mean(mag,2)
    dmag[d,i]=median(sigmag)
   
  endfor
  
endfor

plotout='~/Desktop/cent_stk_sig_v.ps'
plotpdf='~/Desktop/cent_stk_sig_v.pdf'

stop
;window, 0, xsize=600, ysize=500, xpos=1500, ypos=100
;ps_on, plotout, xsize=15, ysize=14
;plot, mean_mag[0,*], dmag[0,*], yrange=[0.002,0.02], color=cgcolor('black'), /nodata, xtitle='Mean Magnitude', $
;  ytitle='Sigma LC (mags)', title='Centaurus Targets - v1', charsize=0.7
;for i=0, ndir-1 do oplot, mean_mag[i,*], dmag[i,*], color=cols1[i], psym=8
;al_legend, ['1stk', '5stk', '10stk'], psym=[8,8,8], color=cols1, /left, textcolor=cgcolor('black')
;ps_off

;spawn, 'convert '+plotout+' '+plotpdf

;stop

ps_on, plotout, xsize=15, ysize=14
plot, vmags[0,*], dmag[0,*], yrange=[0.002,0.02], color=cgcolor('black'), /nodata, xtitle='Vmag', $
  ytitle='Sigma LC (mags)', title='Centaurus Targets - v1', charsize=0.7
for i=0, ndir-1 do oplot, vmags[i,*], dmag[i,*], color=cols1[i], psym=8
al_legend, ['1stk', '5stk', '10stk'], psym=[8,8,8], color=cols1, /left, textcolor=cgcolor('black')
ps_off

spawn, 'convert '+plotout+' '+plotpdf


stop
print, 'end of program'
end