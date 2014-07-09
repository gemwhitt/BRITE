pro ba_ub_temps

Compile_opt idl2

indir='/Users/gemmawhittaker/BRITE/results/ba_ub/ORION/'

filesin=file_search(indir+'*.txt', count=nf)

readcol, filesin[0], jd_ba, temp_ba, format='d,f'

readcol, filesin[1], jd_ub, temp_ub, format='d,f'

plotsym, 0, /fill, 0.4
plot, jd_ub, temp_ub, color=cgcolor('black'), /nodata
oplot, jd_ub, temp_ub, color=cgcolor('purple'), psym=8
oplot, jd_ba, temp_ba, color=cgcolor('green'), psym=8

stop
end