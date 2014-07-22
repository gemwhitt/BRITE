pro brite_autocor1

; Investigate correlations in the BRITE datasets
; Output: Autocorrelation graph

Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.3

sat='BA'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/' ; p1
;indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/HD31237/'  ; p0
filesin=file_search(indir+'*.sav', count=nf)

for f=0, nf-1 do begin

  obj=obj_new('IDL_Savefile', filesin[f])
  obj->restore, 'jd'
  obj->restore, 'data1'
  obj->restore, 'ccd_temp'
  obj->restore, 'vmag'
  obj->restore, 'medimg0'
  
  jd1=jd-jd[0]
  
  nfrm=n_elements(jd)
  
  fname=file_basename(filesin[f],'.sav')
  
  ; calculate the total DN in each image
  totdn=lonarr(nfrm)
  for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
  
  ; calculate the average temperatures in each image
  avgtmp=fltarr(nfrm)
  for img=0, nfrm-1 do avgtmp[img]=average(ccd_temp[*,img])
  
  ; make shorter subsets of the data
  avgtmp1=avgtmp[10000:10199]
  totdn1=totdn[10000:10199]
  n1=n_elements(totdn1)
  
  plotsym, 0, /fill, 0.8
  plot, avgtmp1, totdn1, color=cgcolor('black'), psym=8, /ynozero
  stop
  
  ; DO CROSS-CORRELATIONS....
  print, correlate(avgtmp1, totdn1)
  lag=indgen(200)
  
  result=c_correlate(avgtmp1, totdn1, lag)
  
  plot, result, color=cgcolor('black'), psym=8
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ; DO AUTO-CORRELATIONS.......
  ac=fltarr(n_elements(lag))  
  for i=0, n_elements(lag)-1 do begin 
    nitt=n1-lag[i]
    ac1=fltarr(nitt)
    
    for k=0, nitt-1 do ac1[k]=(totdn1[k] - mean(totdn1)) * (totdn1[k+lag[i]] - mean(totdn1))
        
    ac[i]=total(ac1)*(1./(n1-lag[i]))
  endfor
  
  result=ac/ac[0]  
  
  plot, result, color=cgcolor('black'), psym=8
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;Now what??
  
  stop
  
endfor

end