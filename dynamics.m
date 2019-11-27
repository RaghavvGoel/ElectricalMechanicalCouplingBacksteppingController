function [q_dot_dot,I_dot]  = dynamics(Tau_d,q_dot,q, I, V)
%#codegen
p1=3.473;
p2=0.196;
p3=0.242;
fd1=5.3;
fd2=1.1;

c2=cos(q(2));
s2=sin(q(2));
% alpha=1;

M=[p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];% Inertia matrix
Vm=[-p3*s2*q_dot(2) -p3*s2*(q_dot(1)+q_dot(2));p3*s2*q_dot(1) 0];%Centripetal coriolis matrix
fd=[fd1 0;0 fd2]; % Friction matrix
q_dot_dot=-M\(Vm*q_dot+fd*q_dot-Tau_d);   %EL dynamics

%as L = I2, inv(L) = I2;
I_dot = V - q_dot - I;