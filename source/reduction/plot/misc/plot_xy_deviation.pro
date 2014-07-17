pro plot_xy_deviation

; program to measure the deviation in x and y for each target over the duration of the observations
; get min, max and mean - plot
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.5, /fill

indir='~/BRITE/data/UB/p4/CENTAURUS/1stk/lc_txt/'

outdir='~/BRITE/data/UB/p4/CENTAURUS/1stk/plots/stats/'

filesin=file_search(indir+'*.txt', count=nf)

fileout=outdir+'xy_deviations.txt'
plotout=outdir+'xy_deviations.pdf'

for i=0, nf-1 do begin
  
  readcol, filesin[i], skipline=14, jd, xpsf, ypsf, format='d,x,x,x,x,f,f,x,x'
  
  readcol, filesin[i], skipline=2, name, format='a', numline=1
  
  readcol, filesin[i], skipline=6, vmag, format='a', numline=1
  
  npt=n_elements(jd)
  
  jd1=jd-jd[0]
  jd2=jd1[1:npt-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  gap=[0,gap,npt-1]
  
  max_dz=fltarr(ngap+1)
  med_dz=fltarr(ngap+1)
  sig_dz=fltarr(ngap+1)
  
  for j=0, ngap do begin
    
    iloc=indgen(gap[j+1]+1)
    
    n=n_elements(iloc)
    
    x2=xpsf[1:n-1]
    dx=abs(x2-xpsf[iloc])
    
    y2=ypsf[1:n-1]
    dy=abs(y2-ypsf[iloc])
    
    dz=sqrt(dx^2 + dy^2)
    
    max_dz[j]=max(dz)
    med_dz[j]=median(dz)
    sig_dz[j]=robust_sigma(dz/median(dz))
    
  endfor
  
;window, 0, xsize=600, ysize=500, xpos=1500, ypos=100
;plot, max_dz, color=cgcolor('black'), psym=8

openw, lun, fileout, /get_lun, /append
if i eq 0 then printf, lun, 'max_shift', 'mean_shift', 'sig_shift', format='(a10,x,a10,x,a10)'
printf, lun, min(max_dz), median(med_dz), median(sig_dz), format='(d10.2, x, d10.2, x, d10.2)'
free_lun, lun

up1=med_dz*sig_dz+med_dz
low1=med_dz-med_dz*sig_dz

;window, 1, xsize=600, ysize=500, xpos=2300, ypos=100
tempplot1=outdir+'tempplot.ps'
tempplot2=outdir+'tempplot.pdf'
tempplot3=outdir+'tempplot3.pdf'

ps_on, tempplot1, xsize=15, ysize=14
plot, indgen(ngap+1), med_dz, color=cgcolor('black'), psym=8, yrange=[0,5], xtitle='Orbit', $
  ytitle='Average shift in PSF center', title=name+', '+vmag, charsize=0.8
errplot, indgen(ngap+1), low1, up1, color=cgcolor('black'), width=0.01
ps_off

spawn, 'convert '+tempplot1+' '+tempplot2

if i eq 0 then begin
  
spawn,  'mv '+tempplot2+' '+plotout 

endif else begin
    
  spawn, 'mv '+plotout+' '+tempplot3
  
  spawn, 'PDFconcat -o '+plotout+' '+tempplot2+' '+tempplot3
  
endelse

endfor

print, 'end of program'
end