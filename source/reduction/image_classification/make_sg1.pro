pro make_sg1

; make sg1 files and plots for further analysis
; 
; 
Compile_opt idl2
!p.background=cgcolor('white')

indir='~/BRITE/data/UB/p2/ORION/'

outdir='~/BRITE/data/UB/p2/ORION/sg1/'

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin
  
  hdname=file_basename(filesin[f], '_p2.sav')
  
  obj=obj_new('IDL_Savefile', filesin[f])
  obj->restore, 'jd'
  obj->restore, 'data1'
  obj->restore, 'vmag'
  obj->restore, 'ccd_temp'
  obj->restore, 'iflag'
  obj->restore, 'medimg0'
  
  
  good=where(iflag eq 2, ngood, complement=bad)
  
  nbad=n_elements(bad)
  
  nfrm=ngood
  
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  data1=data1[*,*,good]
  jd=jd[good]
  ccd_temp=ccd_temp[*,good]
  ccd_temp2=fltarr(nfrm)
;  stop
  for jj=0, nfrm-1 do ccd_temp2[jj]=avg(ccd_temp[*,jj])
  ccd_temp=ccd_temp2
  
  jd1=jd-jd[0]
  
  npix=[]
  flx=[]
  
  for j=0, nfrm-1 do begin
  
    data2=reform(data1[*,*,j])
    
    totdn=total(data2)
    
    mostdn=0.5*totdn
    
    data3=reform(data2, xdim*ydim)
    sort1=reverse(sort(data3))
    data3=data3[sort1]
    subtot=fltarr(xdim*ydim)
    for jj=0, xdim*ydim-1 do subtot[jj]=total(data3[0:jj])
    xx=(where(subtot ge mostdn))[0]
    
    ; save xx = number of pixels containing most of the flux + save what most flux is!
    npix=[npix,xx]
    flx=[flx,mostdn]
    
  endfor
  
    plot, npix, flx, color=cgcolor('black'), psym=2, xtitle='Number of pixels', ytitle='50% total flux', charsize=1.0, $
    title=hdname
    
  ; SG1 - images with good PSFs
  sub1=where(npix le 50, nsub1)
  med1=median(flx[sub1])
  sig1=stddev(flx[sub1])
    
  keep1=where(flx[sub1] ge med1-5*sig1 AND flx[sub1] le med1+5*sig1, nkeep1, complement=rej1)
  
  oplot,  npix[sub1[keep1]], flx[sub1[keep1]], color=cgcolor('purple'), psym=2
  
  oplot,  npix[sub1[rej1]], flx[sub1[rej1]], color=cgcolor('pale green'), psym=2
 
 ; SG2 - images with good blank frames
 sub2=where(npix gt 70, nsub2)
 med2=median(flx[sub2])
 sig2=stddev(flx[sub2])
 
 keep2=where(flx[sub2] lt med2+2*sig2, nkeep2, complement=rej2)
 
 oplot,  npix[sub2[keep2]], flx[sub2[keep2]], color=cgcolor('sky blue'), psym=2
 
  oplot,  npix[sub2[rej2]], flx[sub2[rej2]], color=cgcolor('red'), psym=2

  
  print, vmag
  
  stop
  continue
  
  sg1=sub1[keep1]
  nsg1=nkeep1
  
  ; make plot to show which images are sg1
  plotout=outdir+'plots/'+hdname+'_sg1.ps'
  pdfout=outdir+'plots/'+hdname+'_sg1.pdf'
  
  ps_on, plotout, xsize=17, ysize=16
  
  plot, jd1, ccd_temp, color=cgcolor('black'), psym=2, xtitle='Time (days)', ytitle='Temperature', $
    title=hdname, charsize=0.7
  oplot, jd1[sg1], ccd_temp[sg1], color=cgcolor('purple'), psym=2
  
  ps_off

  spawn, 'convert '+plotout+' '+pdfout
  spawn, 'rm '+plotout
  
  ; make save file
  jd=jd[sg1]
  jd1=jd1[sg1]
  data1=data1[*,*,sg1]
  ccd_temp=ccd_temp[sg1]
 
  fileout=outdir+'sav/'+hdname+'_sg1.sav'
  save, filename=fileout, vmag, hdname, jd, data1, ccd_temp
  
  
  ; save stat results to stats file
  pctar=fix(float(nsub1)/float(nfrm)*100.)
  pcgood=fix(float(nkeep1)/float(nfrm)*100.)
  
  openw, lun, outdir+'sg1_stats.dat', /get_lun, /append
  printf, lun, hdname, vmag, pctar, pcgood, med1, sig1, format='(a7,x,f7.2,x,i4,x,i4,x,f7.1,x,f7.1)'
  free_lun, lun
  
endfor


print, 'end of program'
end