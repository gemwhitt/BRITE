pro thrm_tst_results1

; plot residuals from the thermal tests to show whether there is a difference 
; pre-post > 0 means that there is a lower dark current in the annealed image ---> it has worked
; pre-post < 0 means that the annealing made the dark current worse!
; pre-post = 0 means no difference
; 
Compile_opt idl2
;
; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

; input directories
indir='~/BRITE/data/anneal_test_data/sav_files/'

outdir='~/BRITE/results/annealing/date_8nov13/'
outfile=outdir+'plots/pixel_dist.ps'
!p.multi=0
ps_on, outfile, xsize=15, ysize=15

savfiles=file_search(indir+'*.sav', count=nsav)

if nsav eq 0 then stop

atemp=[0.06,1.,10.]
cols=[cgcolor('blue'), cgcolor('dark green'), cgcolor('orange')]

; compare the distrubution of pixel values
for ii=0, 2 do begin

  dates=strarr(nsav)
  temps=fltarr(nsav)
  for i=0, nsav-1 do begin
  
  restore, savfiles[i]
  dates[i]=strmid(fdate,4)
  temps[i]=temp[ii]
  
  if i eq 0 then cghistoplot, data[*,*,ii], binsize=5, backcolorname=cgcolor('white'), axiscolorname=cgcolor('black'), polycolor=cols[i], xrange=[-10,200], $
    charsize=0.7, ytitle='Number of pixels', xtitle='Pixel data value', title='Readout time = '+strmid(strtrim(atemp[ii],2),0,4)+'s' , $
    datacolorname=cols[i], /line_fill, orientation=45 else $
    cghistoplot, data[*,*,ii], binsize=5, polycolor=cols[i], $
    charsize=0.7, /oplot, /line_fill, orientation=45, datacolorname=cols[i]
  endfor
  
  al_legend, dates+' - '+strmid(strtrim(temps,2),0,4), textcolor=cols, /right, charsize=0.6


endfor

ps_off

spawn, 'open '+outfile+' &'
stop
end
