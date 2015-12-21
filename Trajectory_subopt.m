
function [a,b] = Trajectory_subopt(q0,t0,qf,tf,o1)

% This program is for static obstacles

%a=[a0 a1 a2 a3 a4 a5 a6]
%b=[b0 b1 b2 b3 b4 b5 b6]
%xd(t)=a0+a1t+a2t^2+a3t^3+a4t^4+a5t^5+a6t^6.
%yd(t)=b0+b1t+b2t^2+b3t^3+b4t^4+b5t^5+b6t^6.
%t0: initial time
%q0: vehicle initial configuration, q0=[x0 y0 theta0 phi0 v0 dv0]
%tf: ending time
%qf: vehicle ending configuration, qf=[xf yf thetaf phif vf dvf]
%(x,y): position.
%theta: heading angle (deg).
%phi: steering angle (deg).
%v: velocity.
%dv: acceleration.

%vehicle parameters:
%L: the distance between the front and rear axles.

%Obstacles:
%o1:(x1,y1,r1), the obstacle's center and radius. If no obstacle exists, then set r1=0.

format long;

x0=q0(1); y0=q0(2); theta0=q0(3)*pi/180; phi0=q0(4)*pi/180; v0=q0(5); dv0=q0(6);
xf=qf(1); yf=qf(2); thetaf=qf(3)*pi/180; phif=qf(4)*pi/180; vf=qf(5); dvf=qf(6);
x1=o1(1); y1=o1(2); r1=o1(3); L=2; yt=30; yb=10; r0=5;

h1=x0;
h2=v0*cos(theta0);
h3=dv0*cos(theta0)-v0^2*tan(phi0)*sin(theta0)/L;
h4=xf;
h5=vf*cos(thetaf);
h6=dvf*cos(thetaf)-vf^2*tan(phif)*sin(thetaf)/L;

h7=y0;
h8=v0*sin(theta0);
h9=dv0*sin(theta0)+v0^2*tan(phi0)*cos(theta0)/L;
h10=yf;
h11=vf*sin(thetaf);
h12=dvf*sin(thetaf)+vf^2*tan(phif)*cos(thetaf)/L;

L=[1 t0 t0^2 t0^3 t0^4 t0^5;
   0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;
   0 0 2 6*t0 12*t0^2 20*t0^3;
   1 tf tf^2 tf^3 tf^4 tf^5;
   0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;
   0 0 2 6*tf 12*tf^2 20*tf^3];

H1=[h1;h2;h3;h4;h5;h6];
H2=[h7;h8;h9;h10;h11;h12];
M=[t0^6;6*t0^5;30*t0^4;tf^6;6*tf^5;30*tf^4];

syms t ft f1 f2 f3 xp yp;
ft=[1 t t^2 t^3 t^4 t^5];
f1=t^6-ft*inv(L)*M;
f2=ft*inv(L)*H1;
f3=ft*inv(L)*H2;
kx=(xf-x0)/(tf-t0);
ky=(yf-y0)/(tf-t0);
xp=kx*(t-t0)+x0;
yp=ky*(t-t0)+y0;
f1
f2
f3
p1=double(vpa(int(f1^2,t0,tf)));
p2=double(vpa(int(f1*(f2-xp),t0,tf)));
p3=double(vpa(int(f1*(f3-yp),t0,tf)));
p4=double(vpa(int((f2-xp)^2+(f3-yp)^2,t0,tf)));

if(r1<0.01) %no obstacles

a6_star=-p2/p1;
b6_star=-p3/p1;

a=inv(L)*(H1-M*a6_star);
b=inv(L)*(H2-M*b6_star);

a=[a;a6_star]
b=[b;b6_star]

end

%*************************%
%   Collision avoidance   %
%*************************%

if(r1>0.01)
    
    JSubStar_case1=inf;
    JSubStar_case2=inf;
    JSubStar_case3=inf;
    JSubStar_case4=inf;
    
    b6max_case1=-inf;b6min_case1=inf;
    a6max_case2=-inf;a6min_case2=inf;
    a6max_case3=-inf;a6min_case3=inf;
    a6max_case4=-inf;a6min_case4=inf;

    T=0.01; %calculation step.
    t=t0:T:tf;
    k=length(t);

    for i=21:k-20
        ft=[1 t(i) t(i)^2 t(i)^3 t(i)^4 t(i)^5];
        f1=t(i)^6-ft*inv(L)*M;
        f2=ft*inv(L)*H1;
        f3=ft*inv(L)*H2;
        
        %case 1: a6=0    
        
        if(abs(r0+r1)>abs(f2-x1))
            b6max_tmp=sqrt((r0+r1)^2-(f2-x1)^2)/abs(f1)-(f3-y1)/f1;
            b6min_tmp=-sqrt((r0+r1)^2-(f2-x1)^2)/abs(f1)-(f3-y1)/f1;
            if(b6max_case1<b6max_tmp)
                b6max_case1=b6max_tmp;
            end
            if(b6min_case1>b6min_tmp)
                b6min_case1=b6min_tmp;
            end
        end
        
        %case 2: b6=0 
        
        if(abs(r0+r1)>abs(f3-y1))
            a6max_tmp=sqrt((r0+r1)^2-(f3-y1)^2)/abs(f1)-(f2-x1)/f1;
            a6min_tmp=-sqrt((r0+r1)^2-(f3-y1)^2)/abs(f1)-(f2-x1)/f1;
            if(a6max_case2<a6max_tmp)
                a6max_case2=a6max_tmp;
            end
            if(a6min_case2>a6min_tmp)
                a6min_case2=a6min_tmp;
            end
        end   
        
        %case 3: a6=b6
        
        Delta=(f2+f3-x1-y1)^2-2*((f2-x1)^2+(f3-y1)^2-(r0+r1)^2);
        if(Delta>0)
            a6max_tmp=-(f2+f3-x1-y1)/(2*f1)+sqrt(Delta)/(2*abs(f1));
            a6min_tmp=-(f2+f3-x1-y1)/(2*f1)-sqrt(Delta)/(2*abs(f1));
            if(a6max_case3<a6max_tmp)
                a6max_case3=a6max_tmp;
            end
            if(a6min_case3>a6min_tmp)
                a6min_case3=a6min_tmp;
            end
        end  
        
        %case 4: a6=-b6
        
        Delta=(f2-f3-x1+y1)^2-2*((f2-x1)^2+(f3-y1)^2-(r0+r1)^2);
        if(Delta>0)
            a6max_tmp=-(f2-f3-x1+y1)/(2*f1)+sqrt(Delta)/(2*abs(f1));
            a6min_tmp=-(f2-f3-x1+y1)/(2*f1)-sqrt(Delta)/(2*abs(f1));
            if(a6max_case4<a6max_tmp)
                a6max_case4=a6max_tmp;
            end
            if(a6min_case4>a6min_tmp)
                a6min_case4=a6min_tmp;
            end
        end
    end

%
%   For case 1: a6=0
%
    b6Star_case1=-p3/p1;
    a6_case1=0;
    b6_case1=b6Star_case1;
    if(b6Star_case1>b6min_case1&&b6Star_case1<b6max_case1)
        a6_case1=0;
        b6_case1=b6max_case1;
        a=inv(L)*(H1-M*a6_case1);
        b=inv(L)*(H2-M*b6_case1);
        a=[a;a6_case1];
        b=[b;b6_case1];
        xd=a(1)+a(2)*t+a(3)*t.^2+a(4)*t.^3+a(5)*t.^4+a(6)*t.^5+a(7)*t.^6;
        yd=b(1)+b(2)*t+b(3)*t.^2+b(4)*t.^3+b(5)*t.^4+b(6)*t.^5+b(7)*t.^6;
        if((yd<yt)&(yd>yb))
            a6_case1=0;
            b6_case1=b6max_case1;
        else
            a6_case1=0;
            b6_case1=b6min_case1;        
        end
    end
    JSubStar_case1=p1*a6_case1^2+2*p2*a6_case1+p1*b6_case1^2+2*p3*b6_case1+p4;

%
%   For case 2: b6=0
%  
    a6Star_case2=-p2/p1;
    a6_case2=a6Star_case2;
    b6_case2=0;
    if(a6Star_case2>a6min_case2&&a6Star_case2<a6max_case2)
        a6_case2=a6max_case2;
        b6_case2=0;
        a=inv(L)*(H1-M*a6_case2);
        b=inv(L)*(H2-M*b6_case2);
        a=[a;a6_case2];
        b=[b;b6_case2];
        xd=a(1)+a(2)*t+a(3)*t.^2+a(4)*t.^3+a(5)*t.^4+a(6)*t.^5+a(7)*t.^6;
        yd=b(1)+b(2)*t+b(3)*t.^2+b(4)*t.^3+b(5)*t.^4+b(6)*t.^5+b(7)*t.^6;
        if((yd<yt)&(yd>yb))
            a6_case2=a6max_case2;
            b6_case2=0;       
        else
            a6_case2=a6min_case2;
            b6_case2=0;       
        end 
    end
    JSubStar_case2=p1*a6_case2^2+2*p2*a6_case2+p1*b6_case2^2+2*p3*b6_case2+p4;

%
%   For case 3: a6=b6
%    
    a6Star_case3=-(p2+p3)/(2*p1);
    a6_case3=a6Star_case3;
    b6_case3=a6_case3;
    if(a6Star_case3>a6min_case3&&a6Star_case3<a6max_case3)
        a6_case3=a6max_case3;
        b6_case3=a6max_case3;
        a=inv(L)*(H1-M*a6_case3);
        b=inv(L)*(H2-M*b6_case3);
        a=[a;a6_case3];
        b=[b;b6_case3];
        xd=a(1)+a(2)*t+a(3)*t.^2+a(4)*t.^3+a(5)*t.^4+a(6)*t.^5+a(7)*t.^6;
        yd=b(1)+b(2)*t+b(3)*t.^2+b(4)*t.^3+b(5)*t.^4+b(6)*t.^5+b(7)*t.^6;
        if((yd<yt)&(yd>yb))
            a6_case3=a6max_case3;
            b6_case3=a6max_case3;      
        else
            a6_case3=a6min_case3;
            b6_case3=a6min_case3;       
        end  
    end
    JSubStar_case3=p1*a6_case3^2+2*p2*a6_case3+p1*b6_case3^2+2*p3*b6_case3+p4;

%
%   For case 4: a6=-b6
%
    a6Star_case4=-(p2-p3)/(2*p1);
    a6_case4=a6Star_case4;
    b6_case4=-a6_case4;
    if(a6Star_case4>a6min_case4&&a6Star_case4<a6max_case4)
        a6_case4=a6max_case4;
        b6_case4=-a6_case4;
        a=inv(L)*(H1-M*a6_case4);
        b=inv(L)*(H2-M*b6_case4);
        a=[a;a6_case4];
        b=[b;b6_case4];
        xd=a(1)+a(2)*t+a(3)*t.^2+a(4)*t.^3+a(5)*t.^4+a(6)*t.^5+a(7)*t.^6;
        yd=b(1)+b(2)*t+b(3)*t.^2+b(4)*t.^3+b(5)*t.^4+b(6)*t.^5+b(7)*t.^6;
        if((yd<yt)&(yd>yb))
            a6_case4=a6max_case4;
            b6_case4=-a6_case4;          
        else
            a6_case4=a6min_case4;
            b6_case4=-a6_case4;          
        end
    end
    JSubStar_case4=p1*a6_case4^2+2*p2*a6_case4+p1*b6_case4^2+2*p3*b6_case4+p4;
    
    tmp=[JSubStar_case1;JSubStar_case2;JSubStar_case3;JSubStar_case4];
    
    [JSubStar,Idx]=min(tmp)
    
    if(Idx==1)
        a6_star=a6_case1;
        b6_star=b6_case1;
    elseif(Idx==2)
        a6_star=a6_case2;
        b6_star=b6_case2;
    elseif(Idx==3)
        a6_star=a6_case3;
        b6_star=b6_case3;
    elseif(Idx==4)
        a6_star=a6_case4;
        b6_star=b6_case4;        
    end
    
    a=inv(L)*(H1-M*a6_star);
    b=inv(L)*(H2-M*b6_star);
    a=double([a;a6_star]);
    b=double([b;b6_star]);
 end
