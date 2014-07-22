pro brite_fft_noise

; Investigate FFT to determine frequency of systematic noise?
; Output: Power spectra
  
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
   
  im=data1[*,*,0]
      
  plot_image, bytscl(im, 20, 200)
  
  f=fft(im, -1)
  
  pow=abs(f)^2
  
  ;m=indgen(n1)-(float(n1)/2. - 1)
  ;freq=m/
  stop

  plot_image, bytscl(pow, 20, 200)
  
  
  
  stop
   

endfor
  
end