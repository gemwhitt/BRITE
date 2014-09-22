pro move_testdata

; move files (without editing) from original directory to a different different directory
; append filenames with name of satellite...
  
Compile_opt idl2
  
; for selecting files to go into testsets
sat=['UB','BA','LEM','TOR']
nsat=n_elements(sat)

field='CENTAURUS'
  
target=['HD127973','HD129056']  ; for centaurus
;target=['HD35411','HD37128']  ; for orion
  
; which level of file
level='p2'

outdir='~/BRITE/TESTSETS/werner4lc/'+level+'/

for s=0, nsat-1 do begin    ; begin loop over number of satellites

  filein=file_search('~/BRITE/'+sat[s]+'/'+field+'/data/'+level+'/'+target+'*', count=nf)
  
  fname=file_basename(filein, level+'.sav')  
    
  for f=0, nf-1 do spawn, 'cp '+filein[f]+' '+outdir+fname[f]+sat[s]+'.sav'
    
endfor
  
print, 'End of program'
print, 'Analyse testsets!'
  
end

