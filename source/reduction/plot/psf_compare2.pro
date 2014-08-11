pro psf_compare2

sat=['UB','BA','LEM','TOR']
field='CENTAURUS'

; outdir
indir='~/BRITE/'+sat+'/'+field+'/plots/p1/PSFs/'

filesin=file_search(indir+'*t19.ps', count=nf)

fname=file_basename(filesin, '.ps')

outdir='~/BRITE/results/PSF_compare/CENTAURUS/'

; get target id
dashpos=strsplit(fname, '_')  ; this will be a structure e.g.[0,2] is first target
target=strarr(nf)
for i=0, nf-1 do target[i]=strmid(fname[i], 0, dashpos[i,1]-1)

; get unique targets
utar=target[uniq(target, sort(target))]
nutar=n_elements(utar)

; get the temps
temps=strarr(nf)
for i=0, nf-1 do temps[i]=strmid(fname[i], dashpos[i,4]+1)

for i=0, nutar-1 do begin  ; concatenate 4 images of same target from different satellites

  fout=outdir+utar[i]+'_'+temps[0]+'.ps'
  
  ; get the files for this target
  xx=where(target eq utar[i], nxx)
  
  ;for j=0, nxx-1 do begin
  ;  lenname=strlen(filesin[xx[j]])
  ;  newname=strmid(filesin[xx[j]], 0, lenname-2)
  ;  ;spawn, 'convert '+filesin[xx[j]]+' -resize 150% '+newname+'png'
  ;  spawn, 'convert '+filesin[xx[j]]+' '+newname+'ps'
  ;  spawn, 'convert '+newname+'ps'+' '+newname+'png'
  ;  ;stop
  ;endfor
  
  if nxx eq 1 then $
    spawn, 'montage -geometry 275x275>+1+1 '+filesin[xx[0]]+' '+fout
  if nxx eq 2 then $
    spawn, 'montage -geometry 275x275>+1+1 '+filesin[xx[0]]+' '+filesin[xx[1]]+' '+fout
   if nxx eq 3 then $
    spawn, 'montage -tile 2x2 -geometry 275x275>+1+1 '+filesin[xx[0]]+' '+filesin[xx[1]]+' '+filesin[xx[2]]+' '+fout
  if nxx eq 4 then $
    spawn, 'montage -geometry 275x275>+1+1 '+filesin[xx[0]]+' '+filesin[xx[1]]+' '+filesin[xx[2]]+' '+filesin[xx[3]]+' '+fout
  
endfor


print, 'End of program'
end