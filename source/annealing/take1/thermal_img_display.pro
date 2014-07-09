pro thermal_img_display

; read in and compare images using plot_image with appropriate scaling - use this to compare (by-eye) the images before and after annealing

Compile_opt idl2
  
; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')
  
; input directory
indir='~/BRITE/data/anneal_test_data/sav_files/'
  
savfiles=file_search(indir+'*.sav', count=nsav)

count=0
for i=0, nsav-1 do begin
  
  restore, savfiles[i]  ;fdate, temp, data, readtime
  
  print, file_basename(savfiles[i],'.sav')
  
  for j=0, 2 do begin
    window, j, xsize=600, ysize=450, xpos=1500+(j*700), ypos=200
    plot_image, bytscl(data[*,*,j], 50, 5000), color=cgcolor('black'), title=strtrim([j],2)+', temp='+strtrim(temp,2)
  endfor
  stop
  
endfor

stop

end