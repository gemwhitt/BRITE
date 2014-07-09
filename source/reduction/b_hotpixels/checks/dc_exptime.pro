pro dc_exptime   

; program to track increase/decrease in the image medians with exposure time and CCD temp.
;
; 26 Feb 2014
;
; CREATED BECAUSE HP-MAPS DON'T WORK SO WELL FOR LONGER EXPOSURE IMAGES AS THEY DO FOR 1S, 0.1S, ...
; THEREFORE NEED TO UNDERSTAND HOW PIXEL DN CHANGES WITH EXP_TIME, FOR A GIVEN CCD_TEMP
;
Compile_opt idl2
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

indir='~/BRITE/data/UB/p1/'

filesin=file_search(indir+'Orion-CF1-2*.sav', count=nsav)

outdir=



