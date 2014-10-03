pro warsaw_compare

Compile_opt idl2


target=['HD127973_UB.txt','HD127973_TO.txt'] ;,'HD129056']

indir='~/BRITE/ALL_SATS/'

;outdir='~/BRITE/'+sat+'/'+field+'/plots/aper_lc_final/'

filein=file_search(indir+target+'*', count=nf)

if nf eq 0 then stop

plotsym, 0, /fill, 0.7
cols1=[cgcolor('purple'),cgcolor('green'),cgcolor('red'),cgcolor('blue')]

for ff=0, nf-1 do begin
  
  readcol, filein[ff], time, frame, flux, xcom, $
    ycom, temperature, vmag, bmag, resid, $
    format='(d,i,f,f,f,f,f,f,f)'
    
    time1=time-2456836.D
    
    if ff eq 0 then plot, time1, flux, xrange=[0,6], color=cgcolor('black'), psym=8, yrange=[50000,72000], $
      xtitle='Time - 2456836D', ytitle='Flux (DN)', charsize=0.9
    
    oplot, time1, flux, color=cols1[ff+2], psym=8
    
    ; now calculate the error
    xx=where(time1 le 6, nxx)
    print, nxx
    time1=time1[xx]
    flux=flux[xx]
    time2=time1[1:(nxx-1)]
    tdiff=time2-time1
    gap=where(tdiff gt 0.015, ngap)
    gap=[-1,gap,nxx-1]
    
    err1=[]
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      ni=n_elements(iloc)
      
      flux=flux/median(flux)
      mags=-2.5*alog10(flux)
      ;mags=mags/median(mags)
      
      err1=[err1,(stddev(mags[iloc])/sqrt(float(ni)))^2.]
      
    endfor
    
    xx=where(finite(err1) eq 0, nxx, complement=keep)
    err1=err1[keep]
    nn=ngap+1-nxx
    
    err2=sqrt(total(err1, /nan)/float(nn))
    
    print, filein[ff], vmag[0], err2*1000., ' mmag'
    
endfor

plotsym, 0, /fill, 1.
al_legend, ['UB','TOR'], color=cols1[2:3], psym=8, /bottom, /right, textcolor=cgcolor('black')






end