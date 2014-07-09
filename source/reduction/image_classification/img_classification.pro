pro img_classification

Compile_opt idl2

sat='BA'
field='ORION'

indir='~/BRITE/data/'+sat+'/p2/'+field+'/'

outdir='~/BRITE/data/'+sat+'/p2/'+field+'/class/'

statsfile=outdir+'stats/class_stats1.txt

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin

  hdname=file_basename(filesin[f], '_p2.sav')
  
  restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
                      ;medcol1, medimg0, ndead, nsat, simbad_radec, vmag, bmag, $
                      ;parlax, otype, sptype, iflag
                      
  nfrm=n_elements(jd)
  ccd_temp2=fltarr(nfrm)
  for jj=0, nfrm-1 do ccd_temp2[jj]=avg(ccd_temp[*,jj])
   
  good=where(iflag eq 2, ngood, complement=bad)  
  nbad=n_elements(bad)
  
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
 
  jd1=jd-jd[0]
  
  npix=[]
  flx=[]
  
  for j=0, ngood-1 do begin
  
    data2=reform(data1[*,*,good[j]])
    
    totdn=total(data2)
    
    mostdn=0.9*totdn
    
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
  
 ; plot, npix, flx, color=cgcolor('black'), psym=2
  
  ; divvy up:
  ; SG1 - images with good PSFs
  sub1=where(npix le median(npix)+(5*robust_sigma(npix)), nsub1, complement=sub2)
  med1=median(flx[sub1])
  sig1=stddev(flx[sub1])
  
  keep1=where(flx[sub1] ge med1-5*sig1 AND flx[sub1] le med1+5*sig1, nkeep1, complement=rej1)
  
  sg1=[sub1[keep1]]
  sg1r=[sub1[rej1]]
  
  ; SG2 - images with good blank frames
  ;sub2=where(npix gt 70, nsub2)
  med2=median(flx[sub2])
  sig2=robust_sigma(flx[sub2])
  
  keep2=where(flx[sub2] lt med2+1*sig2, nkeep2, complement=rej2)
  
  sg2=[sub2[keep2]]
  sg2r=[sub2[rej2]]
  
  ;sub3=where(npix gt 50 AND npix le 70, nsub3)
  
  ; check if any of the in-betweens have flx ~ median flux - then add these to sub2[keep2]
 ; keep3=where(flx[sub3] ge med1-5*sig1 AND flx[sub3] le med1+5*sig1, nkeep3, complement=rej3)
  
 ; if nkeep3 gt 0 then sg1=[sg1,[sub3[keep3]]]
  
  ; check if any of the in-betweens have flx ~ median flux - then add these to sub2[keep2]
 ; keep4=where(flx[sub3] lt med2+1*sig2, nkeep4, complement=rej4)
  
 ;if nkeep4 gt 0 then sg2=[sg2, [sub3[keep4]]]
 
 ; get last subgroup - neither target - nor background = unknown
  ;if nkeep3 gt 0 then sub3[keep3]=-9
  ;if nkeep4 gt 0 then sub3[keep4]=-9
  ;x3=where(sub3 gt 0, nx3)
  ;if nx3 gt 0 then sg3=sub3[x3]

  ; for PLOTTING..... ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;wset, 0
  plotsym, 0, /fill,0.5
  plotout=outdir+'plots/'+hdname+'_class.ps'
  ps_on, plotout, xsize=17, ysize=28
  !p.multi=[0,1,2,0,0]
  plot, npix, flx, color=cgcolor('black'), xtitle='Number of pixels', ytitle='90% total flux', charsize=0.7, $
    title=hdname, psym=8
      
  oplot,  npix[sg1], flx[sg1], color=cgcolor('purple'), psym=8
  
  oplot,  npix[sg1r], flx[sg1r], color=cgcolor('pale green'), psym=8
  
  oplot,  npix[sg2], flx[sg2], color=cgcolor('sky blue'), psym=8
  
  oplot,  npix[sg2r], flx[sg2r], color=cgcolor('red'), psym=8
  
 ; oplot,  npix[sg3], flx[sg3], color=cgcolor('green'), psym=8
  
  al_legend, ['target', 'blank', 'bad', 'unknown'], textcolor=cgcolor(['purple', 'sky blue', 'red', 'pale green']), /left
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;wset, 1
  plot, jd, ccd_temp2, color=cgcolor('black'), charsize=0.7, xtitle='JD', ytitle='CCD Temp', psym=8
 ; oplot, jd[bad], ccd_temp2[bad], color=cgcolor('red'), psym=8
  oplot, jd[good[sg2r]], ccd_temp2[good[sg2r]], color=cgcolor('red'), psym=8
  oplot, jd[good[sg1]], ccd_temp2[good[sg1]], color=cgcolor('purple'), psym=8
  oplot, jd[good[sg2]], ccd_temp2[good[sg2]], color=cgcolor('sky blue'), psym=8
  oplot, jd[good[sg1r]], ccd_temp2[good[sg1r]], color=cgcolor('green'), psym=8
 ; oplot, jd[good[sg3]], ccd_temp2[good[sg3]], color=cgcolor('green'), psym=8
  
  ps_off
  !p.multi=0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;save iflag - 2 is target, 1 is blank, 0 is bad, 3 is unknown
  iflag[good[sg1r]]=0
  iflag[good[sg2]]=1
  iflag[good[sg2r]]=0
;  iflag[good[sg3]]=3
  
  ; write out stats
  trgt=where(iflag eq 2, ntar)
  blnk=where(iflag eq 1, nblnk)
  bad=where(iflag eq 0, nbad)
;  qmrk=where(iflag eq 3, nq)
  
  tot=n_elements(jd)
  
  pct=float(ntar)/float(tot)*100.
  pcbl=float(nblnk)/float(tot)*100.
  pcba=float(nbad)/float(tot)*100.
;  pcq=float(nq)/float(tot)*100.
  
  openw, lun, statsfile, /get_lun, /append
  printf, lun, 'hdname', 't_frms', 'target', 'blank', 'bad', 'unknown', format='(a7,x,a7,x,a7,x,a7,x,a7,x,a7)'
  printf, lun, hdname, strtrim(tot,2), pct, pcbl, pcba, format='(a7,x,a7,x,f7.3,x,f7.3,x,f7.3)'
  free_lun, lun
 
  ; re-save p2 file
  savout=outdir+'sav/'+file_basename(filesin[f])
  save, filename=savout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
                          medcol1, medimg0, ndead, nsat, simbad_radec, vmag, bmag, $
                          parlax, otype, sptype, iflag, ccd_temp2
  
   
endfor

spawn, 'open '+statsfile+' &'

print, 'end of program'
end