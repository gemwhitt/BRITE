pro combine_pdf

Compile_opt idl2

; Program to combine images into one file with 6 images on a page to show centroid tracking results

indir='~/BRITE/results/UB/pointing/images/'
outdir='~/BRITE/results/UB/pointing/pdf_images/'

dates='20130717'

imagefiles=file_search(indir+'*_'+dates+'.ps', count=ntot)

; get ROIs from image names

fnames=file_basename(imagefiles, '_'+dates+'.ps')
dashpos=strsplit(fnames, '_')
rois=strarr(ntot)
for i=0, ntot-1 do rois[i]=strmid(fnames[i], 0, dashpos[i,2]-1)
uniq_rois=rois[uniq(rois, sort(rois))]

for i=0, n_elements(uniq_rois)-1 do begin
  match, rois, uniq_rois[i], suba, subb, count=count1
  
  ;spawn, '"/System/Library/Automator/Combine PS Pages.action/Contents/Resources/join.py" -o ~/Desktop/FILE.pdf '+imagefiles[suba]
  spawn, 'PDFconcat -o ~/Desktop/outfile.pdf ~/BRITE/results/UB/pointing/images/*.pdf'

  stop
endfor

stop
print, 'End of Program'
end



