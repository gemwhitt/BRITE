pro img_classification2

Compile_opt idl2
  
sat='BA'
field='ORION'
  
indir='~/BRITE/data/'+sat+'/p2/'+field+'/'
  
outdir='~/BRITE/data/'+sat+'/p2/'+field+'/class/'
  
statsfile=outdir+'stats/class_stats.txt
  
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
    
  frac=fltarr(ngood)
  flag=intarr(ngood)
  cumdn=lonarr(ngood)
    
  for j=0, ngood-1 do begin
    
    data2=reform(data1[*,*,good[j]])
    
    data3=reform(data2, xdim*ydim)  ; 1D array of DN
    order1=reverse(sort(data3))     ; order from highest to lowest value of DN
    cen=order1[0]              ; 1D coords of brightest/central pixel
    
    ; add up DN in the top 50 pixels
    cumdn[j]=total(data3[order1[0:49]])
    
    continue
    
    ; make a 5x5 cutout around the central pixel
    xcen=(array_indices(data2, cen))[0]
    ycen=(array_indices(data2, cen))[1]
    if xcen-2 lt 0 then x1=0 else x1=xcen-2
     if xcen+2 gt xdim-1 then x2=xdim-1 else x2=xcen+2
      if ycen-2 lt 0 then y1=0 else y1=ycen-2
       if ycen+2 gt ydim-1 then y2=ydim-1 else y2=ycen+2
       
    
    cutout=data2[x1:x2,y1:y2]
    
    totcut=total(cutout)
    
    avgdn=totcut/25.
    stop
    if avgcut gt 100 then flag[j]=1 else flag[j]=2
    
    data4=data3[order1]
    totd4=total(data4[0:24])
    
    frac[j]=totcut/totd4

  endfor
    
   wset, 0
   plot, cumdn, color=cgcolor('black'), psym=2
   
   med=median(cumdn)
   sig=robust_sigma(cumdn)
   up=4*sig+med
   low=med-(4*sig)
   
   oplot, [0,8000], [med,med], thick=3, color=cgcolor('purple')
   oplot, [0,8000], [low,low], thick=3, color=cgcolor('purple')
   oplot, [0,8000], [up, up], thick=3, color=cgcolor('purple')
   

   stop
   continue
   wset, 1
   plot, flag, color=cgcolor('black'), psym=2
   
   x1=where(flag eq 1, nx1, complement=x2)
   
   for kk=0, n_elements(x2)-1 do begin
    plot_image, bytscl(data1[*,*,good[x2[kk]]], 20, 200)
    
    wait, 0.6
   endfor
   
   stop
   
   
   
   
   
   
   
    ;save iflag - 2 is target, 1 is blank, 0 is bad, 3 is unknown
    iflag[good[sg1r]]=3
    iflag[good[sg2]]=1
    iflag[good[sg2r]]=0
    iflag[good[sg3]]=3
    
    ; write out stats
    trgt=where(iflag eq 2, ntar)
    blnk=where(iflag eq 1, nblnk)
    bad=where(iflag eq 0, nbad)
    qmrk=where(iflag eq 3, nq)
    
    tot=n_elements(jd)
    
    pct=float(ntar)/float(tot)*100.
    pcbl=float(nblnk)/float(tot)*100.
    pcba=float(nbad)/float(tot)*100.
    pcq=float(nq)/float(tot)*100.
    
    openw, lun, statsfile, /get_lun, /append
    printf, lun, 'hdname', 't_frms', 'target', 'blank', 'bad', 'unknown', format='(a7,x,a7,x,a7,x,a7,x,a7,x,a7)'
    printf, lun, hdname, strtrim(tot,2), pct, pcbl, pcba, pcq, format='(a7,x,a7,x,f7.3,x,f7.3,x,f7.3,x,f7.3)'
    free_lun, lun
    
    ; re-save p2 file
    savout=outdir+'sav/'+file_basename(filesin[f])
    save, filename=savout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
      medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
      parlax, otype, sptype, iflag, ccd_temp2
      
      
  endfor
  
  spawn, 'open '+statsfile+' &'
  
  print, 'end of program'
end