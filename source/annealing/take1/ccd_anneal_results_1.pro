pro ccd_anneal_results_1

Compile_opt idl2

; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

; 15th October - first program to analyse new data from thermal testing 
; Compare new results (for annealed CCD, with old ones (pre-anneal) 
; Questions: Is it worth annealing the UB CCD for the purpose of healing HPs - without further degrading the pixels
; 
; Need a global and individual analysis of pixels
; 
; 1. What is the distribution of HPs pre and post annealing 
; 2. What percentage of pixels are HPs/CPs??
; 3. What is the behaviour of HPs are different readout times?
; 4. Compare certain HPs pre and post annealing - start here
; 
; input directory
indir1='~/BRITE/data/anneal_test_data/2013Sept24_ThermalTest_Triumf_Proton_test_CCD/'
indir2='~/BRITE/data/anneal_test_data/2013Oct11_AnnealingTest_Triumf_Proton_test_CCD/'

; save results to....
outdir='~/BRITE/results/annealing_151013/'

;times=['60ms1', '60ms2', '1s1', '1s2', '10s']
times=['1s1']
ntimes=n_elements(times)

;read_temp=['0','10','20','30']
read_temp=['20','30']
ntemp=n_elements(read_temp)

for i=0, ntemp-1 do begin
  
   pre_fits=file_search(indir1+read_temp[i]+'_'+times+'.fits')
   pre_data=mrdfits(pre_fits, 0, header1)
   
   window, 0, xsize=600, ysize=500, xpos=1500, ypos=150
   plot_image, bytscl(pre_data, 50, 150), color=cgcolor('black'), title='Pre Anneal '+strtrim(read_temp[i],2)
   
   post_fits=file_search(indir2+'*_'+times+'.fits')
   post_data=mrdfits(post_fits, 0, header2)
   
   window, 1, xsize=600, ysize=500, xpos=2300, ypos=150
   plot_image, bytscl(post_data, 50, 150), color=cgcolor('black'), title='Post Anneal'

stop 
endfor

print, 'End of Program'
end