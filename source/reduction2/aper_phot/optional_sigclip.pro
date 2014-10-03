pro optional_sigclip

; Do a final sigma clip on the corrected light curves
; 
Compile_opt idl2

sat='UB'

field='CENTAURUS'

target=['HD127973','HD129056']

indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/'

outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/'
chk=file_search(outdir, count=nchk)
if nchk eq 0 then spawn, 'mkdir -p '+outdir

filein=file_search(indir+target+'*.txt', count=nf)

plotsym, 0, /fill, 0.6

for ff=0, nf-1 do begin

  readcol, filein[ff], time, frame, flux, xcom, $
    ycom, temperature, vmag, bmag, resid, $
    format='(d,i,f,f,f,f,f,f,f)'
    
    caldat, time, mon, day, yr
    xx=where(yr lt 2013, nxx, complement=keep)
    
    time=time[keep]
    frame=frame[keep]
    flux=flux[keep]
    xcom=xcom[keep]
    ycom=ycom[keep]
    temperature=temperature[keep]
    resid=resid[keep]
    
  nimg=n_elements(time)
    
  time1=time-time[0]
  time2=time1[1:nimg-1]
  diff=time2-time1
  gap=where(diff gt 0.015, ngap)
  xx=indgen(gap[0])
  gap=[-1,gap,nimg-1]
  cadence=median(diff[xx])  ; in days
  totdur=time1[nimg-1]-time1[0]
  
;  wset, 0
;  plot, time1, flux, color=cgcolor('black'), psym=2, /YNOZERO;, xrange=[40,60]
;  stop
;  
  flux2=flux*0
  ; first get the median of the points in each orbit
  med_orb=fltarr(ngap+1)
  med_time=dblarr(ngap+1)
  for orb=0, ngap do begin
    iloc=indgen(gap[orb+1]-gap[orb])+(gap[orb]+1)
    
    med_orb[orb]=median(flux[iloc])
    
    med_time[orb]=median(time1[iloc])
    
  endfor
  
 ; oplot, med_time, med_orb, color=cgcolor('green'), psym=8
 
;stop
  ; phase fold on the dominate frequency
  ; Define period/freq range
  fmin=1.0/(totdur) > (1.0/15.)  
  
  ;fmin=(1.0/(5.))
  fmax=(1.0/0.01)
  
  scargle, med_time, med_orb, fmin=fmin, fmax=fmax, om, px, signi=signi
  
  f=om/(2.0*!pi)
  
  xx=where(px eq max(px))
  
  peakpow=px[xx]
  peakfreq=f[xx]

  ;plot, time, px, color=cgcolor('black'), psym=2, xrange=[0,f[xx]*1.5]
  
  peakper=1./peakfreq

  ; do alternative phase folding (see BLS)
  ; Phase fold with period
  phase=fltarr(nimg)
  f0=peakfreq[0]
  for j=0, nimg-1 do begin                       ; for each element in the time array, repeat...
    ph=time1[j]*f0                                ; phase=time*freq
    ph=ph-floor(ph)                           ; phase=reciprocal value of phase - i.e. between 0 and 1
    phase[j]=ph                             ; phase calculated with fortran method
  endfor
  
  ;plot, phase, flux, color=cgcolor('black'), psym=8, /ynozero
 ;stop 
  ; smooth over the dominate period by filling in the gaps
  npts2=round(totdur/cadence)
  
  hist=histogram(phase, binsize=0.01, locations=loc, reverse_indices=ri)
  
  fluxsm=flux*0
  
  sig_group=[]
  
  ; calculate the median of the groups in each bin and subtract that value from each point, 
  ; then calculate the sigma of the whole group and clip with median whole group +/- 5*sigma
  for nh=0, n_elements(hist)-1 do begin
    
    if ri[nh+1] eq ri[nh] then continue
   
    iloc=ri[ri[nh]:ri[nh+1]-1]
    
    if n_elements(iloc) eq 0 then continue
    
    med_group=median(flux[iloc])
    
    fluxsm[iloc]=flux[iloc]-med_group
    
    sig_group=[sig_group,sqrt(total((fluxsm[iloc])^2.)/float(n_elements(iloc)))]
      
  endfor
;stop  
 ; plot, phase, fluxsm, color=cgcolor('black'), psym=2
;  
  med_res=median(fluxsm)
  sig_res=median(sig_group)
  
  out=where(fluxsm le med_res-(3*sig_res) OR fluxsm ge med_res+(3*sig_res), nout, complement=keep)
  ;if nout gt 0 then oplot, phase[out], fluxsm[out], psym=8, color=cgcolor('red')
  stop
  
;  plotsym, 0, /fill, 0.8
;  if nout gt 0 then oplot, phase[out], fluxsm[out], color=cgcolor('red'), psym=8
;  stop
;; 
time=time[keep]
flux=flux[keep]
frame=frame[keep]
xcom=xcom[keep]
ycom=ycom[keep]
temperature=temperature[keep]
vmag=vmag[keep]
bmag=bmag[keep]
resid=resid[keep]

 
 plot, time, flux, color=cgcolor('black'), psym=8, /ynozero
 
 
 
 ; write out new file with correction applied
fileout=outdir+file_basename(filein[ff])
openw, lun, fileout, /get_lun
for ii=0, n_elements(keep)-1 do  printf, lun, time[ii], frame[ii], flux[ii], xcom[ii], $
    ycom[ii], temperature[ii], vmag[ii], bmag[ii], resid[ii], $
    format='(d14.6,x,i,x,f,x,f,x,f,x,f7.3,x,f7.3,x,f7.3,x,d9.2)'
free_lun, lun
 
endfor
  
print, 'End of Progam'
print, 'Now plot LCs and do stats'
end
