function [A,b] = lagrangeDerive(n)
% This function derives the equations of motion of an n-pendulum using
% Lagrange equations and a minimal set of coordinates,e.g. thetas 1,2 & 3 for a
% triple pendulum

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
% construct er and etheta vectors for each bar
% er is fixed to the bar and parallel to it
for i = 1:n
    er(1:3,i) = [cos(t(i)); sin(t(i)); 0];
    et(1:3,i) = [-sin(t(i)); cos(t(i)); 0];
end


% velocity of CGs
vG = sym('vG',[3,n],'real'); %initialize space
for j=1:n
    vOprime = sym(zeros(3,1));
    for k=1:j-1 % contributions from bars above jth bar
        vOprime = vOprime + td(k)*l(k)*et(1:3,k);
    end
    % contribution from jth bar
    vG(:,j) = vOprime + td(j)*d(j)*et(1:3,j);
end

% kinetic energy
KE=sym(zeros(n,1));
for i = 1:n
        KE(i) = .5*m(i)*vG(:,i)'*vG(:,i) + .5*Ig(i)*td(i)^2;
end
E_k = sum(KE);

% potential energy
PE = sym(zeros(n,1));
for i=1:n % for each bar
        h=0;
        for j=1:i-1 % for each bar above ith bar, get the height contribution
            h = h + l(j)*cos(t(j)); % contribution from bars above bar of interest
        end
        h=h + d(i)*cos(t(i)); % contribution from bar of interest
        PE(i) = -m(i)*g*h; % for the ith bar
end

E_p = sum(PE);

% the Lagrangian
L = E_k - E_p;

% find n lagrange equations
% diff wrt each theta (for each theta...
% mechanics of it:
% compute dL/dq-dot and then do chain rule d()/dq-dot * dq-dot/dt to
% differentiate wrt time for the first lagrange term
% sources: I adapted the following beautiful lines from Ruina's demo in class
EoM = jacobian(jacobian(L,td)',[t td])*[td tdd]' - jacobian(L, t')';


[A,b] = equationsToMatrix(EoM, tdd);
