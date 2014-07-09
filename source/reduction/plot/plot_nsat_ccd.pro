pro plot_nsat_ccd

; program to read the .dat file produced by median_column.pro and plot nsat versus CCD location
; 
Compile_opt idl2

filein='~/Desktop/orion_sat_imgs.dat'

readcol, filein, file, nsat, xc, yc, format='a,i,i,i'

sort1=sort(xc)


hd=strmid(file, 8, 7)
uhd=hd[uniq(hd, sort(hd))]
nhd=n_elements(uhd)

nsat1=intarr(nhd)

for i=0, nhd-1 do begin
  
  iloc=where(hd eq uhd[i])
  
  nsat1[i]=total(nsat[iloc])
  
endfor

xc=xc[0:14]
yc=yc[0:14]

sort1=sort(xc)

yc=yc[sort1]
nsat1=nsat1[sort1]
xc=xc[sort1]    

xccd=4000
yccd=2700

plotsym, 0, /fill, 1.0
plot, xc, nsat1, color=cgcolor('black'), xrange=[0,xccd], psym=-8

;bp4=bubbleplot(xc, yc, magnitude=nsat1, color=cgcolor('blue'))



stop
print, 'end of program'
end