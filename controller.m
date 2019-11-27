function [e,r,Tau_d,e_I,V,Ic_dot] = controller(t,q,q_dot,I,Ic) 

c2=cos(q(2));
s2=sin(q(2));
alpha=0.5;
k=0.5;

p1=3.473;
p2=0.196;
p3=0.242;
fd1=5.3;
fd2=1.1;

qd=[5*sin(t);10*cos(t)];
qd_dot=[5*cos(t);-10*sin(t)];
qd_ddot=[-5*sin(t);-10*cos(t)];
e=q-qd;     %Tracking Error
e_dot=q_dot-qd_dot;

r=e_dot+alpha*e;         %Filtered Tracking Error

%----------------------------------------------------------%
M=[p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];% Inertia matrix
Vm=[-p3*s2*q_dot(2) -p3*s2*(q_dot(1)+q_dot(2));p3*s2*q_dot(1) 0];%Centripetal coriolis matrix
fd=[fd1 0;0 fd2]; % Friction matrix

Tau_d =-k*r+Vm*q_dot+fd*q_dot+M*qd_ddot-alpha*M*e_dot-Vm*r-e;         % Controller

c2_dot = -s2*q_dot(2);
s2_dot = c2*q_dot(2);

q_ddot=-M\(Vm*q_dot+fd*q_dot-Tau_d);

qd_ddot=[-5*sin(t);-10*cos(t)];
qd_dddot = [-5*cos(t);+10*sin(t)];
e_ddot=q_ddot-qd_ddot;

r_dot=e_ddot+alpha*e_dot;

Id = Tau_d;
M_dot =[2*p3*c2_dot p3*c2_dot;p3*c2_dot 0];
Vm_dot=[-p3*(s2*q_ddot(2)+s2_dot*q_dot(2)) -p3*(s2*(q_ddot(1)+q_ddot(2))+s2_dot*(q_dot(1)+q_dot(2)));p3*(s2*q_ddot(1)+s2_dot*q_dot(1)) 0];
%Vm_dot=[-p3*(s2_dot*q_dot(2)) -p3*(s2_dot*(q_dot(1)+q_dot(2)));p3*(s2_dot*q_dot(1)) 0];
%Id_dot=-k*(-qd_ddot+alpha*e_dot)+Vm_dot*q_dot+M_dot*qd_ddot+M*qd_dddot-alpha*(-M*qd_ddot-M_dot*e_dot)-Vm_dot*r-Vm*(-qd_ddot + alpha*e_dot) - e_dot;
Id_dot=-k*(r_dot)+Vm*q_ddot+Vm_dot*q_dot+fd*q_ddot+M_dot*qd_ddot+M*qd_dddot-alpha*(M*e_ddot+M_dot*e_dot)-Vm_dot*r-Vm*(r_dot) - e_dot;
e_I = I - Ic;
gamma = 0.01;
Ic_dot = (1/gamma)*(-Ic + Id);

V = q_dot + Ic_dot-(-Id) - r;