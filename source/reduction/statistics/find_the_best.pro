pro find_the_best


  Compile_opt idl2
   
  indir='~/BRITE/data/UB/p4_gem3_stk/stats/'
  
  filesin=file_search(indir+'1s_*.txt', count=nfiles)
  
  outdir=indir
  
  fileout=outdir+'best.txt'
  
  
  
  tempflux=fltarr(15,4)
  tempsig=fltarr(15,4)
  flux=fltarr(15)
  sig=fltarr(15)
  
  
  for i=0, nfiles-1 do begin
  
    readcol, filesin[i], name, mag, nsat, pix, epix, maxdn, emaxdn, meanflux, eflux, $
    format='a,f,i,i,f,f,f,f, f'
    
    hdname=name
    vmags=mag
    tempflux[*,i]=meanflux
    tempsig[*,i]=eflux
    
  endfor
  
  for i=0, 14 do begin
    
    xx=(where(tempsig[i,*] eq min(tempsig[i,*])))[0]
    
    flux[i]=tempflux[i,xx]
    sig[i]=tempsig[i,xx]
    
  endfor
  
  openw, lun, fileout, /get_lun
  for i=0, 14 do printf, lun, hdname[i], vmags[i], flux[i], sig[i], format='(a,x,f,x,f,x,f)'
  free_lun, lun
 
stop 
 
  
  print, 'end of program'
  ;  stop
end


