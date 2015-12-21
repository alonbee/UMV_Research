function DrawCir(fig,x,y,r1)
alpha=0:0.1:2*pi-0.1;
cir_x=x+r1*cos(alpha);
cir_y=y+r1*sin(alpha);
plot(fig,cir_x,cir_y,'r');

