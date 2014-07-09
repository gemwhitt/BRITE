pro plot_psf_info 

; program to plot the results from sav_psf_info - also look at plots made in sr_fitting
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 0.8

nstk=10

indir='/Users/gemmawhittaker/BRITE/data/UB/p4/CENTAURUS/'+strtrim(nstk,2)+'stk/psf/'

filesin=file_search(indir+'*.txt', count=nf)

; define arrays for collecting global results for all targets
npix=[]
fwhm=[]
fl=[]
efl=[]
vmags=[]
bvmags=[]

for i=0, nf-1 do begin
  
  fname=file_basename(filesin[i], '_psf.txt')
  
  ; FILE STRUCTURE.....
  ; 
  ;'FWHM',fwhm[0],fwhm[1], format='(a15,x,f7.3,x,f7.3)'
  ;'NUM_PIX_IN_AP',npix[0],npix[1],npix[2], format='(a15,x,i7,x,i7,x,i7)'
  ;'E_NUM_PIX',e_npix[0],e_npix[1],e_npix[2], format='(a15,x,f7.4,x,f7.4,x,f7.4)'
  ;'AVG_FLUX',fl[0],fl[1],fl[2], format='(a15,x,d14.1,x,d14.1,x,d14.1)'
  ;'E_FLUX',e_fl[0],e_fl[1],e_fl[2], format='(a15,x,f7.4,x,f7.4,x,f7.4)'
  ;'VMAG',vmag, format='(a10,x,f7.3)'
  ;'B-V_MAG',bminv, format='(a15,x,f7.3)'
  ;''
  ;max_dn[i], xy_psf[0,i], xy_psf[1,i], format='(i10, x, f7.2, x, f7.2)'
  
  readcol, filesin[i], skipline=8, max_dn, xpsf, ypsf, format='i,f,f,'
  
  
   
  ;window, 0, xsize=600, ysize=550, xpos=1500, ypos=100
  ;plot, xpsf, max_dn, color=cgcolor('black'), psym=8, xtitle='X-Pixel', ytitle='Max DN', charsize=0.8
  
 ; window, 1, xsize=600, ysize=550, xpos=2300, ypos=100
 ; plot, ypsf, max_dn, color=cgcolor('black'), psym=8, xtitle='Y-Pixel', ytitle='Max DN', charsize=0.8
  
  ; append arrays
  readcol, filesin[i], skipline=1, numline=1, npix0, npix1, npix2, format='x,i,i,i'
  npix=[npix,npix1]
  fwhm=[fwhm,npix0]
  
  readcol, filesin[i], skipline=3, numline=1, fl0, fl1, fl2, format='x,l,l,l'
  fl=[fl,fl1]
  
  readcol, filesin[i], skipline=4, numline=1, efl0, efl1, efl2, format='x,f,f,f'
  efl=[efl,efl1]
  
  readcol, filesin[i], skipline=5, numline=1, vmag, format='x,f'
  vmags=[vmags,vmag]
 ; print, vmag
  
  readcol, filesin[i], skipline=6, numline=1, bminv, format='x,f'
  bvmags=[bvmags,bminv]
  
endfor

; plot results - global

wset, 0
plot, vmags, fl, color=cgcolor('black'), psym=8

plot, vmags, efl, color=cgcolor('black'), psym=8

plot, vmags, fwhm, color=cgcolor('black'), psym=8


stop
print, 'end of program'
end