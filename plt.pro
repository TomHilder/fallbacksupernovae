PRO PLT

data=read_ascii('prot.dat')
r=data.field1[0,0:299]
u=reform(data.field1[1,*],300,1001)

set_plot,'PS'
device,filename='test.ps'
plot,r,u[*,0]
oplot,r,u[*,100]
oplot,r,u[*,100]
oplot,r,u[*,200]
oplot,r,u[*,300]
oplot,r,u[*,400]
oplot,r,u[*,500]
oplot,r,u[*,600]
oplot,r,u[*,700]
oplot,r,u[*,800]
oplot,r,u[*,900]
oplot,r,u[*,1000]
device,/close
$mv test.ps ~/public_html/

END 