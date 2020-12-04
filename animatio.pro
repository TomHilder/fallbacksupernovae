PRO ANIMATIO
plot_data=read_ascii("prot.dat")
steps=size(plot_data.field1[0,*],/N_ELEMENTS)/200L
data=reform(plot_data.field1,5,200,steps)

frames=bytarr(300,300,steps)
set_plot,"X"
window,1,title="ANIMATIO",xsize=300,ysize=300
for i=0,(steps-1) do begin
	plot,data[0,*,i],data[1,*,i]
	frames[0,0,i]=tvrd()
endfor
for i=0,(steps-1) do tv,frames[*,*,i]
END