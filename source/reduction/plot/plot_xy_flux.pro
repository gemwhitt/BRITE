pro plot_xy_flux

; PURPOSE: Plot x and y PSF locations versus flux - to get magntiude of the difference
; 
Compile_opt idl2

sat='BA'
field='ORION'
rbin=1

plotsym, 0, /fill, 0.6

indir='~/BRITE/data/'+sat+'/sos_phott'+strtrim(fix(rbin),2)+'/'+field+'/lc_txt/'

filesin=file_search(indir+'*.txt', count=nf)

if nf eq 0 then stop

for ff=0, nf-1 do begin
  
  readcol, filesin[ff], jd, flux, npix, max_dn, ccd_temp2, xcom, ycom, x1, x2, y1, y2, $
    format='d,d,i,i,d,d,d,dd,d,d', comment='#'
    
  jd0=jd-jd[0]
    
  mags=(-2.5)*alog10(flux)
    
  plot, xcom, mags, color=cgcolor('black'), psym=8, /ynozero
    
    
    stop
  
endfor
stop


print, 'end of program'
end