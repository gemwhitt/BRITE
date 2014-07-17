pro non_stk_list

; produce a list for Bert of non-stacked fits files
; 
Compile_opt idl2

indir='~/BRITE/data/UB/fits_data/rasters/ORION/'  ; location of fits files

outdir='~/Desktop/' 

filesin=file_search(indir+'*.fits', count=nf)

fileout=outdir+'orion_nonstk_fits.txt'

for i=0, nf-1 do begin
  
  data=mrdfits(filesin[i], 0, header, /status, /silent)
  
  fname=file_basename(filesin[i])
  
  exp_ttl=(sxpar(header, 'EXP_TTL'))/1000.    ; in seconds
  
  if exp_ttl eq 1. then begin
    
    openw, lun, fileout, /get_lun, /append
    printf, lun, fname, format='(a)' 
    free_lun, lun

  endif
  
endfor

print, 'end of program'

end

