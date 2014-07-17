pro script_orion1

; script to find target info for a group of targets

Compile_opt idl2

indir='~/BRITE/data/UB/roi_raw_sav/CENTAURUS/'

outdir='~/BRITE/results/simbad/scripts/'

outfile=outdir+'centaurus_targets1.txt'

filesin=file_search(indir+'*_p0.sav', count=nfiles)

; get target names
fnames=file_basename(filesin, '_p0.sav')
; get start of ID (HD....)
hd_pos=strpos(fnames, 'HD')
target_id=strarr(nfiles)
for i=0, nfiles-1 do target_id[i]=strmid(fnames[i], hd_pos[i])

str='query id'
str=replicate(str,nfiles) ; makes an array with the same number of elements as nfiles with the value of 'query id'

data1=strarr(2,nfiles)

data1[0,*]=str
data1[1,*]=target_id


openw,lun, outfile, /get_lun   ; for testing contamination
printf,lun,'echo My Simbad script '
printf,lun,'format object f1 "%IDLIST(A;HD):%COO(d;AD):%FLUXLIST(U,B,V,R;N=F):%OTYPELIST(S):%PLX(V[E]):%SP(S;Q)"'
printf,lun,'set radius 0.5m'
printf,lun,data1
printf,lun,'format display'
free_lun,lun

spawn, 'open -a textedit '+outfile+' &'

print, 'End of program'
print, 'Submit '+outfile+' to SIMBAD'

END




