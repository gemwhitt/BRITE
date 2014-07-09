pro make_color_template

; make a color template
; 
plot_colors=[cgcolor('green yellow'), cgcolor('crimson'), cgcolor('aquamarine'), cgcolor('red'), cgcolor('blue'), cgcolor('forest green'), $
  cgcolor('pink'), cgcolor('plum')]


save, filename='~/IDLWorkspace82/BRITE/resource/plot_colors.sav', plot_colors


print, 'end of program'
end

