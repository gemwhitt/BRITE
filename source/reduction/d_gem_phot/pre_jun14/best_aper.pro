pro best_aper

; program to quickly display the best aperture - used for the LC
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 1.0, /fill

nstk='1'

indir='~/BRITE/data/UB/p4/CENTAURUS/'+nstk+'stk/lc_txt/'

aper=strarr(30)
target=strarr(30)

filesin=file_search(indir+'*.txt', count=nf)
  
for i=0, nf-1 do begin
  
  readcol, filesin[i], skipline=9, numline=1, ap, format='x,a', delimiter=':'
    ap=strtrim(ap,2)
    
  readcol, filesin[i], numline=1, name, format='a'
  
  aper[i]=ap
  target[i]=name
  
endfor

fileout=indir+'best_ap.dat'

openw, lun, fileout, /get_lun
for i=0, nf-1 do printf, lun, target[i], aper[i], format='(a10,x,a10)'
free_lun, lun

spawn, 'open '+fileout+' &' 
    
end