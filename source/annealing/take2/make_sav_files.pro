pro make_sav_files

Compile_opt idl2

; take fits files and make .sav files - take average of 2 files if they have the same times
; need average tempetarures for this....
; 
; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

; input directories
indir=file_search('~/BRITE/data/anneal_test_data/*_Annealing*/', count=ndir)

outdir='/Users/gemmawhittaker/BRITE/data/anneal_test_data/sav_files/'

dirname=file_basename(indir)
;stop
for i=0, ndir-1 do begin
  
  fdate=strmid(dirname[i], 0, 9) ; save this in result file
  
  print, fdate
  
  atimes=['60','1s','10']
  
  data=lonarr(4048,2672,3)
  
  temp=fltarr(3)
  
  for j=0, 2 do begin
    
    print, atimes[j]
    
    fitsfiles=file_search(indir[i]+'/20_'+atimes[j]+'*.fits', count=nfits) ; should be (2x1s, 2x60ms, 1x10s)
    
    ; get times
    fnames=file_basename(fitsfiles, '.fits')   
    ftimes=strmid(fnames, 3)
    
    if nfits eq 0 OR nfits gt 2 then stop
    
    if nfits gt 1 then begin  ; create an average dataset
      
      data1=mrdfits(fitsfiles[0], 0, header, /silent)
      data2=mrdfits(fitsfiles[1], 0, header, /silent)
      
      ; check file dimensions match
      if (size(data1, /dim))[0] ne (size(data2, /dim))[0] OR (size(data1, /dim))[1] ne (size(data2, /dim))[1] then stop
      
      ; make averages
      data[*,*,j]=(float(data1+data2))/2.
     
    endif else data[*,*,j]=mrdfits(fitsfiles[0], 0, header, /silent) 
        
    print, median(data[*,*,j])
    
    ; get temp from .txt file
    tempfile=(file_search(indir[i]+'/Temperatures_20.txt', count=temp_chk))[0]
    
    if temp_chk eq 1 then begin
      readcol, tempfile, time, t1, t2, t3, t4, format='a,f,f,f,f'
      
      time1=strmid(time, 0, 2)
      
      match, time1, atimes[i], suba, subb, count=count1
      
      if count1 eq 0 then stop
      
      all_temps=[t1[suba],t2[suba],t3[suba],t4[suba]]
      
      temp[j]=avg(all_temps)

    endif else temp[j]=20.
    
  endfor
  
  readtime=[0.06,1,10]
  
  short_date=strmid(fdate, 4)
  
  savefile=outdir+short_date+'.sav'
  save, filename=savefile, fdate, temp, data, readtime
  
;  window, 0, xpos=1500, ypos=200
;  plot_image, bytscl(data[*,*,0], 50, 5000)
;  stop
;  plot_image, bytscl(data[*,*,1], 50, 5000)
; stop
;  plot_image, bytscl(data[*,*,2], 50, 5000)
;  stop
;stop  
endfor

print, 'end of program'
end