function DrawTri(fig,L,cl,x,y,theta)
p1x=x+(L+1)*cos(theta);
p1y=y+(L+1)*sin(theta);
p2x=x+(L-1)*cos(theta+2*pi/3);
p2y=y+(L-1)*sin(theta+2*pi/3);
p3x=x+(L-1)*cos(theta-2*pi/3);
p3y=y+(L-1)*sin(theta-2*pi/3);
plot(fig,[p1x p2x],[p1y p2y],cl,[p2x p3x],[p2y p3y],cl,[p3x p1x],[p3y p1y],cl);

