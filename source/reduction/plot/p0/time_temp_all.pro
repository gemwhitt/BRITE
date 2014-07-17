pro time_temp_all

; Plot time versus temperature for the whole observing run
; Highlight where files differ, e.g. if raster sizes change, or similar.
; 
Compile_opt idl2

sat='BA'
field='ORION'

; input directory
indir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/data/raw_sav/'

; input files
filesin=file_search(indir+'*.sav', count=nf)

; output directory
outdir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/plots/'

; output file
fileout=outdir+sat+'_'+field+'_time_temp'

; get filename
fname=strmid(file_basename(filesin, '.sav'), 0, 7)

uname=fname[uniq(fname, sort(fname))]
nu=n_elements(uname)

undefine, filesin
undefine, nf

for ii=0, 0 do begin  ; do only once
  
  filesin=file_search(indir+uname[ii]+'*.sav', count=nf)
  
  time=[]
  temps=[]
  size=[]
  
  for jj=0, nf-1 do begin
    
    obj=obj_new('IDL_Savefile', filesin[jj])
    obj->restore, 'jd'
    obj->restore, 'ccd_temp'
    obj->restore, 'data1'
    
    xdim=(size(data1, /dim))[0]
    
    npt=n_elements(jd)
    
    size=[size, replicate(xdim, npt)]
    
    time=[time,jd]
    
    avg_temp=fltarr(npt)
    
    for kk=0, npt-1 do avg_temp[kk]=average(ccd_temp[*,kk])
    
    temps=[temps,avg_temp]
    
  endfor
  
  npts=n_elements(time)
  
  ; sort in order of jd
  sort1=sort(time)
  
  time=time[sort1]
  temps=temps[sort1]
  size=size[sort1]
  
  ; convert start and end times into actual dates
  caldat, time, mon, day, yr
  dates=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
  
  time=time-2454000D
  
  plotsym, 0, /fill, 0.2
  ps_on, fileout+'.ps', xsize=15, ysize=13
  plot, time, temps, color=cgcolor('black'), /nodata, title=sat+' - '+field+' - '+dates[0]+' to '+dates[npts-1], $
    xtitle='Time (JD - 2454000)', ytitle='Average CCD temperature', charsize=0.7
  oplot, time[where(size eq 32)], temps[where(size eq 32)], psym=8, color=cgcolor('orange red')
  oplot, time[where(size eq 24)], temps[where(size eq 24)], psym=8, color=cgcolor('slate grey')
  al_legend, ['32x32','24x24'], psym=8, colors=[cgcolor('orange red'),cgcolor('slate grey')], /left, textcolor='black', $
    /box, outline_color=cgcolor('black'), charsize=0.7
  ps_off
  
  spawn, 'convert '+fileout+'.ps '+fileout+'.pdf'

  stop
endfor

stop





print, 'End of program'
end