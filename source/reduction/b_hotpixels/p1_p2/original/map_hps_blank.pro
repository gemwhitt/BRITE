pro map_hps_blank

; program to map the HPs using the empty frames - by a given range of index locations
; 
Compile_opt idl2
;; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

; input directory
indir='~/BRITE/data/UB/p1/CENTAURUS/'

outdir='~/BRITE/data/UB/reduction/hot_pixel_maps/CENTAURUS/'

filesin=file_search(indir+'*.sav', count=nsav)
fname=file_basename(filesin, '_p1.sav')

nfrm=63
iloc=indgen(nfrm)+615

if nsav eq 0 then stop

for fi=0, nsav-1 do begin

  restore, filesin[fi]  ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                        ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
                        
       
  data=data1[*,*,iloc]
  
  temp1=ccd_temp[iloc]
 
  xdim=(size(data, /dim))[0]
  ydim=(size(data, /dim))[1]
  
  ; make an average frame
  x=fltarr(xdim, ydim)
  for i=0, nfrm-1 do x=x+data[*,*,i]
  x=x/float(nfrm)
  
  ; cycle through the frames record the pixel value and ccd_temp
  for i=529, xdim*ydim-1 do begin
    
    pix=array_indices(x,i)
    datapix=data[pix[0], pix[1], *]
    
    plot, temp1, datapix, color=cgcolor('black'), psym=2
    
    stop
    
    plotsym, 0, /fill, 0.8
    plot, temp1, datapix, color=cgcolor('black'), psym=8
    
    cghistoplot, temp1, binsize=1, locations=loc, histdata=result, reverse_indices=ri, $
      mininput=floor(min(temp1))
    
    nloc=n_elements(loc)
    pixcor=fltarr(nloc)
    
    for j=0, nloc-1 do begin
      
      values=datapix[ri[ri[j]:ri[j+1]-1]]
      pixcor[j]=median(values)
      
    endfor
    
    if i eq 0 then begin
       
    ffcor=pixcor
    
    endif else begin
      

      ffcor=[[ffcor],[pixcor]]
    endelse
    
    stop
  endfor
  
  
 
endfor


stop

print, 'end of program'
end