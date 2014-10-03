pro concat_lcs, sat, field, target

; combine LCs from same star but different raster sizes
; 
Compile_opt idl2

;sat='BA'
;;
;field='CENTAURUS'
;
;target='HD127963'

indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/'

outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/'

movedir=outdir+'separate'
chk=file_search(movedir, count=nchk)
if nchk eq 0 then spawn, 'mkdir -p '+movedir

filein=file_search(indir+target+'*_*.txt', count=nf)

fname=file_basename(filein, '.txt')
dashpos=strsplit(fname, '_')
fname2=fname
for ff=0, nf-1 do fname2[ff]=strmid(fname[ff], 0, dashpos[ff,1]-1)

ufname=fname2[uniq(fname2, sort(fname2))]
nu=n_elements(ufname)

for u=0, nu-1 do begin
  
  files=file_search(indir+ufname[u]+'*_*.txt', count=nf2)
  
  time1=[]
  frame1=[]
  flux1=[]
  mags1=[]
  bkgd1=[]
  bkgd_err1=[]
  xcom1=[]
  ycom1=[]
  temperature1=[]
  vmag1=[]
  bmag1=[]
  resid1=[]
  medimg1=[]
  
  for f=0, nf2-1 do begin
    readcol, files[f], time, frame, flux, mags, bkgd, bkgd_err, xcom, $
      ycom, temperature, vmag, bmag, resid, medimg, $
      format='(d,i,f,f,f,f,f,f,f,f,f,f,i)'
      
    time1=[time1,time]
    frame1=[frame1,frame]
    flux1=[flux1,flux]
    mags1=[mags1,mags]
    bkgd1=[bkgd1,bkgd]
    bkgd_err1=[bkgd_err1,bkgd_err]
    xcom1=[xcom1,xcom]
    ycom1=[ycom1,ycom]
    temperature1=[temperature1,temperature]
    vmag1=[vmag1,vmag]
    bmag1=[bmag1,bmag]
    resid1=[resid1,resid]
    medimg1=[medimg1,medimg]
    
    spawn, 'mv '+files[f]+' '+movedir   

  endfor
  
  ; check for uniques  
  utime=time1[uniq(time1, sort(time1))]
  
  if n_elements(utime) ne n_elements(time1) then stop
  
  sort1=sort(time1)
  
  time=time1[sort1]
  frame=frame1[sort1]
  flux=flux1[sort1]
  mags=mags1[sort1]
  bkgd=bkgd1[sort1]
  bkgd_err=bkgd_err1[sort1]
  xcom=xcom1[sort1]
  ycom=ycom1[sort1]
  temperature=temperature1[sort1]
  vmag=vmag1[sort1]
  bmag=bmag1[sort1]
  resid=resid1[sort1]
  medimg=medimg1[sort1]
  
  nimg=n_elements(time)
  
  
  fileout=outdir+ufname[u]+'.txt'
  openw, lun, fileout, /get_lun
  for ii=0, nimg-1 do  printf, lun, time[ii], frame[ii], flux[ii], xcom[ii], $
    ycom[ii], temperature[ii], vmag[ii], bmag[ii], resid[ii], $
    format='(d14.6,x,i,x,f,x,f,x,f,x,f7.3,x,f7.3,x,f7.3,x,d9.2)'
  free_lun, lun
    
endfor

end


