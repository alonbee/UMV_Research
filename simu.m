clear;clc;
format long;

%% initial and final conditions
t0=0; q0=[0 15 0 0 10 0];
tf=8; qf=[100 15 0 0 10 0];
o1=[15 15 3];% center (o1,o2)   radius=3
L=2;

%% call planning algorithms
[a,b]=Trajectory_subopt(q0,t0,qf,tf,o1);



%% Data processing
t=t0:0.01:tf;  
xd=a(1)+a(2)*t+a(3)*t.^2+a(4)*t.^3+a(5)*t.^4+a(6)*t.^5+a(7)*t.^6; %x(t)
yd=b(1)+b(2)*t+b(3)*t.^2+b(4)*t.^3+b(5)*t.^4+b(6)*t.^5+b(7)*t.^6; %y(t)

alpha=0:0.1:2*pi-0.1;
ox=o1(1)+o1(3)*cos(alpha);
oy=o1(2)+o1(3)*sin(alpha);

%% figures
figure(1)
hold;
hd(1)=plot(xd,yd,'b');
hd(2)=plot(xd(1),yd(1),'o');
len=length(xd);
hd(3)=plot(xd(len),yd(len),'x');
hd(4)=plot([0 100],[10 10],'m-');
hd(5)=plot([0 100],[20 20],'m--');
hd(6)=plot([0 100],[30 30],'m-');
hd(7)=plot(ox,oy,'r');
grid on;
legend(hd([1 2 3 4 5 7]),'trajectory','q_0','q_f','road edge','lane divider','obstacle');
title('Trajectory');
xlabel('xd(m)');
ylabel('yd(m)');
axis equal
hold

dotxd=a(2)+2*a(3)*t+3*a(4)*t.^2+4*a(5)*t.^3+5*a(6)*t.^4+6*a(7)*t.^5;
dotyd=b(2)+2*b(3)*t+3*b(4)*t.^2+4*b(5)*t.^3+5*b(6)*t.^4+6*b(7)*t.^5;
vd=sqrt(dotxd.^2+dotyd.^2);
figure(2)
plot(t,vd)
grid on;
title('velocity');
ylim([0 40]);   %% 设置y的范�?
xlabel('t(sec)');
ylabel('v(m/sec)');

thetad=atan2(dotyd,dotxd)*180/pi;
figure(3)
plot(t,thetad)
grid on;
title('heading angle');
xlabel('t(sec)');
ylabel('\theta(deg)');

ddotxd=2*a(3)+6*a(4)*t+12*a(5)*t.^2+20*a(6)*t.^3+30*a(7)*t.^4;
ddotyd=2*b(3)+6*b(4)*t+12*b(5)*t.^2+20*b(6)*t.^3+30*b(7)*t.^4;
dotthetad=(ddotyd.*dotxd-dotyd.*ddotxd)./vd.^2;
phid=atan(L*dotthetad./vd)*180/pi;
figure(4)
plot(t,phid)
grid on;
title('steering angle');
xlabel('t(sec)');
ylabel('\phi(deg)');

kt=(ddotyd.*dotxd-dotyd.*ddotxd)./vd.^3;
figure(5)
plot(t,kt)
grid on;
title('curvature');
xlabel('t(sec)');
ylabel('curvature');

%% Make avi animation

fig=figure(6);
set(fig,'DoubleBuffer','on');
ax=gca;
rect = get(fig,'Position'); 
rect(1:2) = [0 0];
mov = avifile('static_obs.avi','compression','None','fps',100);
for k=1:length(t)
    plot(xd(1:k),yd(1:k),'b');
    hold;
    DrawTri(ax,4,'b',xd(k),yd(k),thetad(k));
    DrawCir(ax,o1(1),o1(2),o1(3));
    axis equal;
    grid on;
    hold;
    F = getframe(fig, rect);
    mov = addframe(mov,F);
end
mov = close(mov);