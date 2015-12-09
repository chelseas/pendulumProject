function fourBarLinkage(t0, td0, p)
% fourBarDerive(); % Was used to write function fourBarRHS()
n=3; p.n=n;
% set constant values (m's, l's, etc. to EOM)
p.d = .5*ones(n,1); p.Ig = ones(n,1); p.g=10;
% initial conditions are important
% solve the ODE
options = odeset('relTol',1e-6,'AbsTol',1e-6);
tspan = linspace(0,5,500);
z0 = [t0; td0];
[tout,zout] = ode45(@(t,y)fourBarRHS(t,y,p), tspan, z0, options);
plotStuff(tout,zout,3,p.l);
totalE = energy(tout, zout, p, 'Four Bar Linkage');

function zdot = fourBarRHS(t,z,p)
% unpack
[t1, t2, t3] = assignstuff(z(1:3)); [td1, td2, td3] = assignstuff(z(4:6));
[m1, m2,m3] = assignstuff(p.m); [Ig1, Ig2, Ig3] = assignstuff(p.Ig);
[d1, d2,d3] = assignstuff(p.d); [l1, l2,l3] = assignstuff(p.l);
g=p.g;

A = [ - Ig1 - d1^2*m1*cos(t1)^2 - d1^2*m1*sin(t1)^2 - l1*m3*cos(t1)*(d3*cos(t3) + l1*cos(t1) + l2*cos(t2)) - l1*m3*sin(t1)*(d3*sin(t3) + l1*sin(t1) + l2*sin(t2)) - l1*m2*cos(t1)*(d2*cos(t2) + l1*cos(t1)) - l1*m2*sin(t1)*(d2*sin(t2) + l1*sin(t1)), - Ig2 - l2*m3*cos(t2)*(d3*cos(t3) + l1*cos(t1) + l2*cos(t2)) - l2*m3*sin(t2)*(d3*sin(t3) + l1*sin(t1) + l2*sin(t2)) - d2*m2*cos(t2)*(d2*cos(t2) + l1*cos(t1)) - d2*m2*sin(t2)*(d2*sin(t2) + l1*sin(t1)), - Ig3 - d3*m3*cos(t3)*(d3*cos(t3) + l1*cos(t1) + l2*cos(t2)) - d3*m3*sin(t3)*(d3*sin(t3) + l1*sin(t1) + l2*sin(t2)), - l1*sin(t1) - l2*sin(t2) - l3*sin(t3), l1*cos(t1) + l2*cos(t2) + l3*cos(t3);...
                                                                                                       - l1*m3*cos(t1)*(d3*cos(t3) + l2*cos(t2)) - l1*m3*sin(t1)*(d3*sin(t3) + l2*sin(t2)) - d2*l1*m2*cos(t1)*cos(t2) - d2*l1*m2*sin(t1)*sin(t2),                                                                       - Ig2 - d2^2*m2*cos(t2)^2 - d2^2*m2*sin(t2)^2 - l2*m3*cos(t2)*(d3*cos(t3) + l2*cos(t2)) - l2*m3*sin(t2)*(d3*sin(t3) + l2*sin(t2)),                           - Ig3 - d3*m3*cos(t3)*(d3*cos(t3) + l2*cos(t2)) - d3*m3*sin(t3)*(d3*sin(t3) + l2*sin(t2)),              - l2*sin(t2) - l3*sin(t3),              l2*cos(t2) + l3*cos(t3);...
                                                                                                                                                                                           - d3*l1*m3*cos(t1)*cos(t3) - d3*l1*m3*sin(t1)*sin(t3),                                                                                                                                                   - d3*l2*m3*cos(t2)*cos(t3) - d3*l2*m3*sin(t2)*sin(t3),                                                                       - Ig3 - d3^2*m3*cos(t3)^2 - d3^2*m3*sin(t3)^2,                            -l3*sin(t3),                           l3*cos(t3);...
                                                                                                                                                                                                                                     -l1*sin(t1),                                                                                                                                                                                             -l2*sin(t2),                                                                                                         -l3*sin(t3),                                      0,                                    0;...
                                                                                                                                                                                                                                      l1*cos(t1),                                                                                                                                                                                              l2*cos(t2),                                                                                                          l3*cos(t3),                                      0,                                    0];
 

b =   [g*m3*(d3*sin(t3) + l1*sin(t1) + l2*sin(t2)) + m2*(l1*cos(t1)*td1^2 + d2*cos(t2)*td2^2)*(d2*sin(t2) + l1*sin(t1)) - m2*(l1*sin(t1)*td1^2 + d2*sin(t2)*td2^2)*(d2*cos(t2) + l1*cos(t1)) + g*m2*(d2*sin(t2) + l1*sin(t1)) - m3*(d3*cos(t3) + l1*cos(t1) + l2*cos(t2))*(l1*sin(t1)*td1^2 + l2*sin(t2)*td2^2 + d3*sin(t3)*td3^2) + m3*(d3*sin(t3) + l1*sin(t1) + l2*sin(t2))*(l1*cos(t1)*td1^2 + l2*cos(t2)*td2^2 + d3*cos(t3)*td3^2) + d1*g*m1*sin(t1);...
                                                                                                       m3*(d3*sin(t3) + l2*sin(t2))*(l1*cos(t1)*td1^2 + l2*cos(t2)*td2^2 + d3*cos(t3)*td3^2) - m3*(d3*cos(t3) + l2*cos(t2))*(l1*sin(t1)*td1^2 + l2*sin(t2)*td2^2 + d3*sin(t3)*td3^2) + g*m3*(d3*sin(t3) + l2*sin(t2)) + d2*m2*sin(t2)*(l1*cos(t1)*td1^2 + d2*cos(t2)*td2^2) - d2*m2*cos(t2)*(l1*sin(t1)*td1^2 + d2*sin(t2)*td2^2) + d2*g*m2*sin(t2);...
                                                                                                                                                                                                                                                                                  d3*g*m3*sin(t3) + d3*m3*sin(t3)*(l1*cos(t1)*td1^2 + l2*cos(t2)*td2^2 + d3*cos(t3)*td3^2) - d3*m3*cos(t3)*(l1*sin(t1)*td1^2 + l2*sin(t2)*td2^2 + d3*sin(t3)*td3^2);...
                                                                                                                                                                                                                                                                                                                                                                                             l1*cos(t1)*td1^2 + l2*cos(t2)*td2^2 + l3*cos(t3)*td3^2;...
                                                                                                                                                                                                                                                                                                                                                                                             l1*sin(t1)*td1^2 + l2*sin(t2)*td2^2 + l3*sin(t3)*td3^2];
 
 

 uks = A\b;
 % state: t1 t2 t3 td1 td2 td3
 zdot = [z(4:6); uks(1:3)];

function plotStuff(tout,zout,n,l)
t = zout(:,1:n); % thetas
% construct r vectors again
rL = zeros(3,n,length(tout));
er = zeros(3,n);
for k = 1:length(tout)
    for i = 1:n
        er(1:3,i) = [cos(t(k,i)); sin(t(k,i)); 0];
    end
    for i = 1:n
        for j = i:n
            rL(1:3,j,k) = rL(1:3,j,k) + l(i)*er(:,i);
        end
    end
end
len = sum(l);
figure(); hold on; axis([-len len -len len]);
title('Four Bar Linkage Motion');
for k = 1:length(tout)
    x = 0;
    y = 0;
    for l = 1:n
        x = [x, rL(1,l,k)];
        y = [y, rL(2,l,k)];
    end
    h = plot(y,-x,'x-k');
    pause(.01);
    delete(h);
end
%plot final positions of bars
x = 0;
y = 0;
for l = 1:n
    x = [x, rL(1,l,end)];
    y = [y, rL(2,l,end)];
end
plot(y,-x,'x-k');
% plot trajectories of ends of 1st and second bars
plot(squeeze(rL(2,1,:)), squeeze(-rL(1,1,:)),'b');
plot(squeeze(rL(2,2,:)), squeeze(-rL(1,2,:)),'m');
plot(squeeze(rL(2,3,:)), squeeze(-rL(1,n,:)),'g');
 
function fourBarDerive()
% NOTE: copied from newton_pendulum_derive
% this funciton outputs matrices that can be solved at each timestep for
% theta1dd, theta2dd, and theta3dd, Rx and Ry
%
n=3; % hard coded for 4 bars (really 3)
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
t = sym('t',[n,1],'real');
td = sym('td',[n,1],'real');
tdd = sym('tdd',[1,n],'real');
syms Ry Rx real; % constraint forces where 3rd bar is pinned
% construct constraint force vector
R = Rx*ei + Ry*ej;
% construct er and etheta vectors for each bar
% er is fixed to the bar and parallel to it
for i = 1:n
    er(1:3,i) = [cos(t(i)); sin(t(i)); 0];
    et(1:3,i) = [-sin(t(i)); cos(t(i)); 0];
end
% construct position vectors to ends of rods and to centers of gravity of rods
% Note: vector to end of nth rod won't be used
rL = sym(zeros(3,n)); % intialize rL matrix
rG = sym(zeros(3,n)); % intialize rG matrix
for i = 1:n
    for j = i:n
        rL(1:3,j) = rL(1:3,j) + l(i)*er(:,i);
    end
end

% vectors from origin to cg of each bar
% for the first  bar
rG(1:3,1) = d(1)*er(:,1);
for i=2:n % for each subsequent bar
    rG(1:3,i) = rL(1:3,i-1) + d(i)*er(:,i);
end

% construct acceleration vectors
% contribution from each bar above the bar of interest
a = sym(zeros(3,n));
for i = 1:n-1
    for j = i+1:n
        a(1:3,j) = a(1:3,j) + (l(i)*tdd(i)*et(:,i) - l(i)*(td(i)^2)*er(:,i));
    end
end
% contribution from bar of interest
for i = 1:n
    if i==3
        aSpec = a(1:3,i) + (l(i)*tdd(i)*et(:,i) - l(i)*(td(i)^2)*er(:,i));
    end
    a(1:3,i) = a(1:3,i) + (d(i)*tdd(i)*et(:,i) - d(i)*(td(i)^2)*er(:,i));
end

% n AMB's
for i = 1:n
    % for this angular mom balance
    M=0; % reset moment
    H=0; % reset change in ang mom
    for j = i:n % the ith AMB has link i:n e.g. 2nd has links 2 and 3
        if i>1 % only subtract vectors if we're not at the origin
            rGmod = rG(:,j) - rL(:,i-1); % construct vector from point of angular mom bal to cg we're dealing with ATM
            rLmod = rL(:,j) - rL(:,i-1); % construct vector from point of angular mom bal to bar end we're dealing with ATM
        else
            rGmod = rG(:,j);
            rLmod = rL(:,j);
        end
        % add up moment and dh/dt contributions from each bar, i:n
        M = M + cross(rGmod,m(j)*g*ei);
        if j==3 % add constraint force if we are on bar 3
            M= M + cross(rLmod,R);
        end
        H = H + cross(rGmod,m(j)*a(:,j)) + Ig(j)*tdd(j)*ek; % accel is wrt fixed frame, but r vector is still wrt point of AMB
    end
    AMB(i) = M(3) == H(3);
end


[A,b] = equationsToMatrix([AMB,aSpec(1)==0,aSpec(2)==0],[tdd,Rx,Ry]);
%A = simplify(A); b=simplify(b);
uks = A\b;
