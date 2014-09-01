pro make_testsets

; make test sets for 2 or more satellite datasets - use raw/p1/p2 data

Compile_opt idl2

; for selecting files to go into testsets
;sat=['UB','BA','LEM','TOR']
sat='TOR'
nsat=n_elements(sat)

field='CENTAURUS'

target=['HD127973','HD129056']  ; for centaurus
;target=['HD35411','HD37128']  ; for orion

; which level of file
level='p2'

; observation dates
oday=[1,30]
omon=[6,7]
oyr=[2014,2014]

otemp='*' ; temp range? or all temps = *

for s=0, nsat-1 do begin    ; begin loop over number of satellites
  
    indir='~/BRITE/'+sat[s]+'/'+field+'/data/'+level+'/'
    
    outdir='~/BRITE/TESTSETS/werner4lc/'+level+'/
    ; check outdir exists and if not - make it
    chk=file_search(outdir, count=nchk)
    if nchk eq 0 then spawn, 'mkdir -p '+outdir
    
    filein=file_search(indir+target+'*.sav', count=nfile)
    
    for f=1, nfile-1 do begin
      
      restore, filein[f]  ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
      ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
      
      fname=file_basename(filein[f], '_p2.sav')
      
      ; convert time range to JD
      jds=julday(omon[0],oday[0],oyr[0])
      jde=julday(omon[1],oday[1],oyr[1])
      
      ; select the range of images using dates and temps
      ind1=where(jd ge jds AND jd lt jde, nind1)
      
      stop
      
      if nind1 eq 0 then continue
      
      if otemp ne '*' then stop
      
      ; make new arrays
      jd=jd[ind1]
      ccd_temp=ccd_temp[*,ind1]
      medimg0=medimg0[ind1]
      roi_loc=roi_loc[*,ind1]
      medcol1=medcol1[ind1]
      medcol2=medcol2[ind1]
      data1=data1[*,*,ind1]
      
      plot, jd-jd[0], ccd_temp[0,*], color=cgcolor('black'), psym=2
     stop       
      ; save new arrays in output directory
      fieldname=strmid(field, 0, 3)
      satname=strmid(sat[s], 0, 2)
      fileout=outdir+satname+'_'+fieldname+'_'+fname+'.sav'
      
      save, filename=fileout, jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag
      
    endfor
  
endfor

print, 'End of program'
print, 'Analyse testsets!'

end

