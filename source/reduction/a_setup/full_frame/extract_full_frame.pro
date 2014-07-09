pro extract_full_frame

; created 17/06/13 by GNW
; Read in full frame fits files for each satellite and display image with info
; prepare for data analysis
; 
indir='~/BRITE/data/'
sat='UB'
filedate='20130703'

fitsfiles=file_search(indir+sat+'/full_frame/'+filedate+'*.fits', count=nfiles)

for i=0, nfiles-1 do begin
  
  fits_info, fitsfiles[i], n_ext=next, extname=extname, /silent
  
  data=mrdfits(fitsfiles[i], next[0], header)
  
  ;window, 0
  tvscl, data
       
  stop
  
endfor

stop
end
