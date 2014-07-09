pro combine_months

; program to combine observing months - from p1 files
; 
Compile_opt idl2

indir='~/BRITE/data/UB/p1/ORION/months/'
outdir='~/BRITE/data/UB/p1/ORION/'

filesin=file_search(indir+'*.sav', count=nf)

fname=file_basename(filesin, '_p1.sav')

hdpos=strpos(fname, 'HD')
hdname=strarr(nf)
for i=0, nf-1 do hdname[i]=strmid(fname[i], hdpos[i])

; find unique HDNames
uhdname=hdname[uniq(hdname, sort(hdname))]
nu=n_elements(uhdname)
if nu ne 15 then stop

for i=0, nu-1 do begin
  
  iloc=where(hdname eq uhdname[i], nf2)
    
  if nf2 ne 6 then stop
  
  for j=0, nf2-1 do begin
    
    restore, filesin[iloc[j]] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                              ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype
                              
    nfrm=n_elements(jd)
                              
    ; get actual dates
    caldat, jd, mon, day, yr, hr, min, sec
    
    yr=strmid(strtrim(yr,2), 2, 2)
    mon=strtrim(mon,2)
    for k=0, nfrm-1 do if strlen(mon[k]) eq 1 then mon[k]='0'+mon[k]   
    day=strtrim(day,2)
    for k=0, nfrm-1 do if strlen(day[k]) eq 1 then day[k]='0'+day[k]
    hr=strtrim(hr,2)
    for k=0, nfrm-1 do if strlen(hr[k]) eq 1 then hr[k]='0'+hr[k]
    min=strtrim(min,2)
    for k=0, nfrm-1 do if strlen(min[k]) eq 1 then min[k]='0'+min[k]
    sec=strtrim(sec,2)
    dotpos=strpos(sec, '.')
    for k=0, nfrm-1 do sec[k]=strmid(sec[k], 0, dotpos[k])
    for k=0, nfrm-1 do if strlen(sec[k]) eq 1 then sec[k]='0'+sec[k]
    
    obs_dates=day+mon+yr
    obs_times=hr+min+sec
    
    if j eq 0 then begin
      
      odates=obs_dates
      otimes=obs_times
      roi=roi_name
      enum=exp_num
      radec=ra_dec
      juldat=jd
      data=data1
      rdim=roi_dim
      xcen=xc
      ycen=yc
      ctemp=ccd_temp
      etime=exp_time
      ettl=exp_ttl
      medcol=medcols
      medim0=medimg0
      medim1=medimg1
      nded=ndead
      nsatu=nsat
    
    endif else begin
      
      odates=[odates,obs_dates]
      otimes=[otimes,obs_times]
      roi=[roi,roi_name]
      enum=[enum,exp_num]
      radec=[[radec],[ra_dec]]
      juldat=[juldat,jd]
      data=[[[data]],[[data1]]]
      rdim=[[rdim],[roi_dim]]
      xcen=[xcen,xc]
      ycen=[ycen,yc]
      ctemp=[ctemp,ccd_temp]
      etime=[etime,exp_time]
      ettl=[ettl,exp_ttl]
      medcol=[[medcol],[medcols]]
      medim0=[medim0,medimg0]
      medim1=[medim1,medimg1]
      nded=[nded,ndead]
      nsatu=[nsatu,nsat]
           
    endelse
        
  endfor
  
  roi_name=roi_name[0]
  exp_num=enum
  ra_dec=radec
  jd=juldat
  data1=data
  roi_dim=rdim
  xc=xcen[0]
  yc=ycen[0]
  ccd_temp=ctemp
  exp_time=etime
  exp_ttl=ettl
  medcols=medcol
  medimg0=medim0
  medimg1=medim1
  ndead=nded
  nsat=nsatu

  save, filename=outdir+uhdname[i]+'_p1.sav', roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
    medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype
  
endfor

print, 'end of program'
end