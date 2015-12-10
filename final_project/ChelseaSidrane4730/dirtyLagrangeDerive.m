function [A,b] = dirtyLagrangeDerive(n)
% This function derives the equations of motion of an n-pendulum using
% links with 3 DoF (x,y,theta) and acceleration constraints at the hinge
% joints
% constants
syms g real
%
% set up coordinate system
ei = [1 0 0]'; ej = [0 1 0]'; ek = [0 0 1]';
% d's = distance from hinges to c.g.s
% l's = distance from hinges to ends
% t's = absolute angles to each pendulum from vertical
m = sym('m',[n,1],'real'); % mass
Ig = sym('Ig',[n,1],'real'); % inertia
d = sym('d',[n,1],'real');
l = sym('l',[n,1],'real');
%
t = sym('t',[1,n],'real'); % it is important that these are 1xn for lagrange
td = sym('td',[1,n],'real'); % it is important that these are 1xn for lagrange
tdd = sym('tdd',[1,n],'real'); % it is important that these are 1xn for lagrange
x =  sym('x',[1,n],'real');
xd =  sym('xd',[1,n],'real');
xdd =  sym('xdd',[1,n],'real');
y = sym('y',[1,n],'real');
yd = sym('yd',[1,n],'real');
ydd = sym('ydd',[1,n],'real');
% constraint forces
Rx = sym('Rx',[1,n],'real');
Ry = sym('Ry',[1,n],'real');
q = [x, y, t]; qd = [xd, yd, td];  qdd = [xdd, ydd, tdd];  

% construct er and etheta vectors for each bar
% er is fixed to the bar and parallel to it
for i = 1:n
    er(1:3,i) = [cos(t(i)); sin(t(i)); 0];
    et(1:3,i) = [-sin(t(i)); cos(t(i)); 0];
end

for i = 1:n
    KE(i) = .5*m(i)*(xd(i)^2 + yd(i)^2) + .5*Ig(i)*td(i)^2;
    PE(i) = -m(i)*g*x(i);
end
% the Lagrangian
L = sum(KE) - sum(PE);

% Lagrange equations in one line
EoM = jacobian(jacobian(L,qd)',[q qd])*[qd qdd]' - jacobian(L, q')';

% Dirty Q's
for i=1:n-1
Qx(i) = Rx(i) - Rx(i+1);
Qy(i) = Ry(i) - Ry(i+1);
Qt(i) = -d(i)*et(:,i)'*(Rx(i)*ei + Ry(i)*ej) + (l(i)-d(i))*et(:,i)'*(-Rx(i+1)*ei - Ry(i+1)*ej);
end
Qx(n) = Rx(n);
Qy(n) = Ry(n);
Qt(n) = -d(n)*et(:,n)'*(Rx(n)*ei + Ry(n)*ej);
Q = [Qx'; Qy'; Qt'];

elegant = EoM==Q;

% accel constraints
aG = sym(zeros(3,n));
for i=1:n
    aG(:,i) = xdd(i)*ei + ydd(i)*ej;
end
% t stands for top, b for bottom of bar
for i=1:n
at(:,i) = aG(:,i) + er(:,i)*d(i)*td(i)^2 - d(i)*tdd(i)*et(:,i); 
ab(:,i) = aG(:,i) - er(:,i)*(l(i)-d(i))*td(i)^2 + (l(i)-d(i))*tdd(i)*et(:,i); 
end

inelegant = at(1:2,1) == [0;0]; % fixed top
for i=1:n-1
    inelegant = [inelegant; ab(1:2,i) == at(1:2,i+1)];
end

[A,b] = equationsToMatrix([elegant;inelegant],[xdd,ydd,tdd,Rx,Ry]);

    
    
    
    