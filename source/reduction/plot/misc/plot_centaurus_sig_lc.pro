pro plot_centaurus_sig_lc

; plot sigma light curve for 1stk, 5stk, and 10stk
;
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 1.0, /fill

nstk=['1','5','10']
  
indir='~/BRITE/data/UB/p4/CENTAURUS/'+nstk+'stk/lc_txt/'
ndir=n_elements(indir)
  
outdir='~/Desktop/'

fileout=outdir+'lc_sigmas.txt'
plotout=outdir+'lc_sigma.pdf'

cols1=cgcolor(['blue', 'green', 'red'])

vmags=fltarr(ndir,30)
sigs=fltarr(ndir,30)
nps=fltarr(ndir,30)

for d=0, ndir-1 do begin
  
  filesin=file_search(indir[d]+'*.txt', count=nf)
  
  for i=0, nf-1 do begin
  
    readcol, filesin[i], skipline=6, numline=1, vmag, format='a'
    vmag=strtrim(vmag,2)
    eqpos=strpos(vmag, '=')
    vmag=float(strmid(vmag, eqpos+1))
    
    readcol, filesin[i], skipline=12, numline=1, sig, format='x,a', delimiter=':'
    sig=strtrim(sig,2)
    pcpos=strpos(sig, '%')
    sig=float(strmid(sig, 0, pcpos))
    
    readcol, filesin[i], skipline=10, numline=1, npix, format='x,a', delimiter=':'
    npix=strtrim(npix,2)
    plpos=strpos(npix, '+')
    npix=float(strmid(npix, 0, plpos-1))
    
    vmags[d,i]=vmag
    sigs[d,i]=sig/100.
    nps[d,i]=npix
    
  endfor
  
endfor



window, 0, xsize=600, ysize=500, xpos=100, ypos=100
plot, vmags[0,*], nps[0,*], yrange=[0,100], color=cgcolor('black'), /nodata
for i=0, ndir-1 do oplot, vmags[i,*], nps[i,*], color=cols1[i], psym=8  

stop

sigs[1,*]=sigs[1,*]/sqrt(5)
sigs[2,*]=sigs[2,*]/sqrt(10)


window, 1, xsize=600, ysize=500, xpos=2300, ypos=100
plot, vmags[0,*], sigs[0,*], yrange=[0,0.15], color=cgcolor('black'), /nodata
for i=0, ndir-1 do oplot, vmags[i,*], sigs[i,*], color=cols1[i], psym=8


stop
print, 'end of program'
end