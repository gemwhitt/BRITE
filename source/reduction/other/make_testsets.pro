pro make_testsets

  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  
  Compile_opt idl2

; 10th Feb 2014

Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 0.9

; make test sets using only "good" observations - use this to characterise the data

indir='~/BRITE/data/UB/p2/'

outdir='~/BRITE/data/UB/testing/saturated/savfiles/'

g=1

; number of orbits to keep in test set
xorb=16

cf=['*1-2','*1-7']

targets=['34085','35468','39801']

filesin=file_search(indir+cf[g]+'*'+targets+'*.sav', count=nsav)
print, nsav

for i=0, nsav-1 do begin
  
  fname=file_basename(filesin[i], '_p2.sav')
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
                      ;simbad_mags, parlax, otype, sptype, exp_time, exp_ttl, medcols
                      
                      
                      xx=where(ccd_temp le 26, nxx)
  
  nfrm=n_elements(jd)
  
  expnum=indgen(nfrm)
  
  jd1=jd-jd[0]
  
  keep=where(jd1 ge 0.0 AND expnum ge 340, nkeep)
  
  jd=jd[keep]
  jd1=jd1[keep]
  data1=data1[*,*,keep]
  ccd_temp=ccd_temp[keep]
  if n_elements(medcols) gt 0 then medcols=medcols[*,keep]
  
  ; keep X number of orbits, e.g. xorb=6
  jd2=jd1[1:nkeep-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  gap=[0,gap,nkeep-1]
  
  iloc=indgen(gap[xorb]+1) 
  
  jd=jd[iloc]
  jd1=jd1[iloc]
  data1=data1[*,*,iloc]
  ccd_temp=ccd_temp[iloc]
  if n_elements(medcols) gt 0 then medcols=medcols[*,iloc]
  
  ; get the magnitude
  vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
  
  fileout=outdir+fname+'_p2.sav'

  ;if n_elements(medcols) gt 0 then save, filename=fileout, jd, jd1, data1, ccd_temp, medcols, roi_dim, simbad_mags else $
  save, filename=fileout, jd, jd1, data1, ccd_temp, roi_dim, simbad_mags
  ;stop
  ;save, filename=fileout, jd, jd1, data1, ccd_temp, simbad_mags, p2_trend
  
endfor


print, 'end of program;'
end

