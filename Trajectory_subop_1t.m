function [a,b]=Trajectory_subopt(q0,t0,qf,tf,o1,L)
format long;

H1=[q0(1) q0(5)*cos(q0(3)) q0(6)*cos(q0(3))-(q0(5)^2*tan(q0(4))*sin(q0(3)))/L qf(1) qf(5)*cos(qf(3)) qf(6)*cos(qf(3))-qf(5)^2*tan(qf(4))*sin(qf(3))/L]';
H2=[q0(2) q0(5)*sin(q0(3)) q0(6)*sin(q0(3))+q0(5)^2*tan(q0(4))*sin(q0(3))/L qf(2) qf(5)*sin(qf(3)) qf(6)*sin(qf(3))+qf(5)^2*tan(qf(4))*sin(qf(3))/L]';
M=[t0^6 6*t0^5 30*t0^4 tf^6 6*tf^5 30*tf^4]';
L=[1 t0 t0^2 t0^3 t0^4 t0^5;0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;0 0 2 6*t0 12*t0^2 20*t0^3;1 tf tf^2 tf^3 tf^4 tf^5;0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;0 0 2 6*tf 12*tf^2 20*tf^3];
syms t  %%符号工具  查具体用法
f=[1 t t^2 t^3 t^4 t^5];
f1=t^6-(f/L)*M; 
f2=(f/L)*H1;
f3=(f/L)*H2;

%%static obstacles此处设vy=0 vx=5;
vy=0;
vx=0;
funmax=((sqrt(4*o1(3)^2-(f3-o1(2)-vy*(t-t0))^2)/abs(f1))-(f2-o1(1)-vx*(t-t0))/f1);
funmin=(-(sqrt(4*o1(3)^2-(f3-o1(2)-vy*(t-t0))^2)/abs(f1))-(f2-o1(1)-vx*(t-t0))/f1);


a61=0;
for t1=t0+1:0.01:tf-1
    a=subs(funmax,t1);
  %  plot(a,t1)
   % hold on
    if a>a61
    a61=subs(funmax,t1);
    end
end

a62=0;
for t1=t0+1:0.01:tf-1
    a=subs(funmin,t1);
    if a<a62
    a62=subs(funmin,t1);
    end
end
x1=((qf(1)-q0(1))/(tf-t0))*(t-t0);
y1=((qf(2)-q0(2))/(tf-t0))*(t-t0);
P1=int(f1^2,t0,tf);    %%积分方法不确定  极值没有求
P2=int(f1*(f2-x1),t0,tf);
P3=int(f1*(f2-y1),t0,tf);
P4=int((f2-x1)^2+(f3-y1)^2,t0,tf);
 J1=P1*(a61+P2/P1)^2+P1*(P3/P1)^2+P4-(P2^2+P3^2)/P1;
J2=P1*(a62+P2/P1)^2+P1*(P3/P1)^2+P4-(P2^2+P3^2)/P1;
if J1<J2
   a6=a61;
else
    a6=a62;
end
 

aa=inv(L)*(H1-M*a6);
bb=inv(L)*(H2);
a=[aa;a6];
b=[bb;0];

end

    
    
    


