pro check_num_pts

; Program for observations which are split due to incompatible raster sizes...
; check that ALL .fits files have been saved in one of p0 or p0b.
; 
Compile_opt idl2

fitsdir='/Users/gemmawhittaker/BRITE/BA/ORION/data/raw_fits/rasters/'

indir='/Users/gemmawhittaker/BRITE/BA/ORION/data/raw_sav/'

fitsfiles=file_search(fitsdir+'*.fits', count=nfits)  ; total number of fits files

savfiles=file_search(indir+'*.sav', count=nsav)

; get uniq target names
fname=strmid(file_basename(savfiles, '.sav'), 0, 7)
un_fname=fname[uniq(fname, sort(fname))]
n_un=n_elements(un_fname)



for un=0, n_un-1 do begin
  
  ; check number of files for this target
  tfiles=file_search(indir+un_fname[un]+'*.sav', count=nt)
  
  time=[]
  
  for f=0, nt-1 do begin
    
    obj=obj_new('IDL_Savefile', tfiles[f])
    obj->restore, 'jd'
    
    time=[time,jd]
    
  endfor
  
  ntime=n_elements(time)
  
  ; match ntime to nfits
  if ntime ne nfits then stop 
  
  endfor
  
print, 'End of program - all files OK'
end
