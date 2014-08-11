pro time_temp_all

; Plot time versus temperature for the whole observing run
; Highlight where files differ, e.g. if raster sizes change, or similar.
; 
Compile_opt idl2

sat='UB'
field='ORION'

; input directory
indir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/data/raw_sav/HD31237/'

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
  
  ; find out how many sized rasters
  usize=size[uniq(size, sort(size))]
  nsize=n_elements(usize)
  if nsize gt 2 then stop
  
  ; convert start and end times into actual dates
  caldat, time, mon, day, yr
  dates=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
  
  time=time-2454000D
  
  ; set up some colours for plotting
  cols1=[cgcolor('orange red'),cgcolor('slate grey')]
  
  plotsym, 0, /fill, 0.2
  ps_on, fileout+'.ps', xsize=15, ysize=13
  plot, time, temps, color=cgcolor('black'), /nodata, title=sat+' - '+field+' - '+dates[0]+' to '+dates[npts-1], $
    xtitle='Time (JD - 2454000)', ytitle='Average CCD temperature', charsize=0.7
for rs=0, nsize-1 do oplot, time[where(size eq usize[rs])], temps[where(size eq usize[rs])], psym=8, color=cols1[rs]
  al_legend, [strtrim(usize,2)], psym=8, colors=cols1, /left, textcolor='black', $
    /box, outline_color=cgcolor('black'), charsize=0.7
  ps_off
  
  spawn, 'convert '+fileout+'.ps '+fileout+'.pdf'

  stop
endfor

stop





print, 'End of program'
end