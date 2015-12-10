% Final Project
% Chelsea Sidrane, crs325@cornell.edu
% MAE 4730
% find equations of motion for triple pendulum 3 different ways
% simulate the motion of a triple pendulum and a four bar linkage
% simulate a pendulum of arbitrary length
function FPwrite()
%% Solve the triple pendulum three ways, and animate them on top of each other
n=3; p.n=n; % set number of links
% Initial Conditions
% set constant values (m's, l's, etc. to EOM)
p.m = 3*ones(n,1); p.l = ones(n,1); p.d = .5*ones(n,1); p.Ig = ones(n,1); p.g=10;
% set common tspan n ICs n options stuff
tspan = linspace(0,30,3000);
%t0 = (pi/2)*ones(n,1)+[0; rand(n-3,1); 2*rand(2,1)]; % interesting initial
%conditions for experimenting
vals = linspace(0,sqrt(2*pi),n)';
t0 = vals.^2 + (pi/2)*ones(n,1);
td0 = zeros(n,1);
z0 = [t0; td0];
options = odeset('relTol',1e-9,'AbsTol',1e-9);
%%
% get the equations of motion one way
[A,b] = newton_pendulum_derive(n); %returns all symbolic equations
writeODE(A,b,0); % write RHS file

% solve the ODE
[tout,zout] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
%plotStuff(tout,zout,n,p.l);

%%
% get the equations of motion a second way
[A,b] = lagrangeDerive(n); %returns all symbolic equations
writeODE(A,b,0); % write RHS file

% solve the ODE
[tout2,zout2] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
%plotStuff(tout,zout,n,p.l);

%%
% get the equations of motion a third way
[A,b] = dirtyLagrangeDerive(n); %returns all symbolic equations
writeODE(A,b,1); % write RHS file

%z0 = [z0(1:3)+.00001; z0(4:6)]; % let the initial conditions deviate just a smidge!
% solve the ODE
[tout3,zout3] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
%plotStuff(tout3,zout3,n,p.l);

%%
% compare
plotTripleStuff(tout,zout,zout2,zout3,n,p.l);
pause(3);
totalE = energy(tout, zout, p,'AMB');
pause(4);
totalE = energy(tout2, zout2, p, 'Clean Lagrange');
pause(4);
totalE = energy(tout3, zout3, p, 'Dirty Lagrange');
pause(4);
clear p
%% Limiting Case of Triple Pendulum
n=3; p.n=n; % set number of links
% Initial Conditions
% set constant values (m's, l's, etc. to EOM)
% here the mass is set really large, so it should act like a single rigid
% rod
p.m = 3*ones(n,1) + [zeros(n-1,1); 1000]; p.l = ones(n,1); p.d = .5*ones(n,1); p.Ig = ones(n,1); p.g=10;
% set common tspan n ICs n options stuff
tspan = linspace(0,5,500);
t0 = (pi/2)*ones(n,1); % linspace(0,2,n)';
td0 = zeros(n,1);
z0 = [t0; td0];
options = odeset('relTol',1e-8,'AbsTol',1e-8);
%
% get the equations of motion
[A,b] = lagrangeDerive(n); %returns all symbolic equations
writeODE(A,b,0); % write RHS file
%
% solve the ODE
[tout4,zout4] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
plotStuff(tout4,zout4,n,p.l);
pause(3);
totalE = energy(tout4, zout4, p, 'Limiting Case: Triple Acts Like Single');
pause(4);
clear p
%% Four Bar Linkage
% first with normal IC's
p.m = ones(3,1); p.l = [1;1;1]; 
t0 = [0 + pi/3; pi/2; pi + pi/3]; 
td0 = [0;0;0];
fourBarLinkage(t0, td0, p);
pause(4);
% Next with random IC's
p.m = rand(n,1); p.l = ones(n,1)+rand(n,1); 
t0 = 3.5*rand(n,1)+[.5;1;-3/2];
td0 = zeros(3,1);
fourBarLinkage(t0, td0, p);
pause(4);
clear all;
%% N-Pendulum
n=10; p.n=n; % set number of links
% Initial Conditions
% set constant values (m's, l's, etc. to EOM)
% here the mass is set really large, so it should act like a single rigid
% rod
p.m = ones(n,1); p.l = ones(n,1); p.d = .5*ones(n,1); p.Ig = ones(n,1); p.g=10;
% set common tspan n ICs n options stuff
tspan = linspace(0,10,1000);
t0 = (pi/2)*ones(n,1);
td0 = linspace(0,1,n)';
z0 = [t0; td0];
options = odeset('relTol',1e-8,'AbsTol',1e-8);
%
% get the equations of motion
[A,b] = lagrangeDerive(n); %returns all symbolic equations
writeODE(A,b,0); % write RHS file
%
pause(.1);
% solve the ODE
[tout5,zout5] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
plotStuff(tout5,zout5,n,p.l);
pause(4);
totalE = energy(tout5, zout5, p, 'Clean Lagrange');

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

figure(); hold on; axis([-n n -n n]); axis square;
title('N Pendulum Motion');
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
% plot trajectories of ends of 1st and second and nth bars
plot(squeeze(rL(2,1,:)), squeeze(-rL(1,1,:)),'b');
plot(squeeze(rL(2,2,:)), squeeze(-rL(1,2,:)),'m');
plot(squeeze(rL(2,n,:)), squeeze(-rL(1,n,:)),'g');
% plot final positions of bars
x = 0;
y = 0;
for l = 1:n
    x = [x, rL(1,l,end)];
    y = [y, rL(2,l,end)];
end
plot(y,-x,'x-k');


function plotTripleStuff(tout,zout,zout2,zout3,n,l)
t = zout(:,1:n); % thetas
t2 = zout2(:,1:n); % thetas
t3 = zout3(:,1:n); % thetas
% construct r vectors again
rL = zeros(3,n,length(tout));
er = zeros(3,n);
rL2 = zeros(3,n,length(tout));
er2 = zeros(3,n);
rL3 = zeros(3,n,length(tout));
er3 = zeros(3,n);
%first
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
% second
for k = 1:length(tout)
    for i = 1:n
        er2(1:3,i) = [cos(t2(k,i)); sin(t2(k,i)); 0];
    end
    for i = 1:n
        for j = i:n
            rL2(1:3,j,k) = rL2(1:3,j,k) + l(i)*er2(:,i);
        end
    end
end
% third
for k = 1:length(tout)
    for i = 1:n
        er3(1:3,i) = [cos(t3(k,i)); sin(t3(k,i)); 0];
    end
    for i = 1:n
        for j = i:n
            rL3(1:3,j,k) = rL3(1:3,j,k) + l(i)*er3(:,i);
        end
    end
end

figure(); hold on; axis([-n n -n n]); axis square;
title('N Pendulum Motion');
for k = 1:length(tout)
    x = 0;
    y = 0;
    x2 = 0;
    y2 = 0;
        x3 = 0;
    y3 = 0;
    for l = 1:n
        x = [x, rL(1,l,k)];
        y = [y, rL(2,l,k)];
        x2 = [x2, rL2(1,l,k)];
        y2 = [y2, rL2(2,l,k)];
        x3 = [x3, rL3(1,l,k)];
        y3 = [y3, rL3(2,l,k)];
    end
    h = plot(y,-x,'x-c');
    h2 = plot(y2,-x2,'o-r');
    h3 = plot(y3,-x3,'.-b');
    pause(.02);
    delete(h); delete(h2); delete(h3);
end
% plot final positions of bars
x = 0;
y = 0;
for l = 1:n
    x = [x, rL(1,l,end)];
    y = [y, rL(2,l,end)];
end
plot(y,-x,'x-c');
%
x2 = 0;
y2 = 0;
for l = 1:n
    x2 = [x2, rL2(1,l,end)];
    y2 = [y2, rL2(2,l,end)];
end
plot(y2,-x2,'o-r');
%
x3 = 0;
y3 = 0;
for l = 1:n
    x3 = [x3, rL3(1,l,end)];
    y3 = [y3, rL3(2,l,end)];
end
plot(y3,-x3,'.-b');
% plot trajectories of ends of 1st and second and nth bars
plot(squeeze(rL(2,1,:)), squeeze(-rL(1,1,:)),'r');
plot(squeeze(rL(2,2,:)), squeeze(-rL(1,2,:)),'m');
plot(squeeze(rL(2,n,:)), squeeze(-rL(1,n,:)),'g');
%
plot(squeeze(rL2(2,1,:)), squeeze(-rL2(1,1,:)),'r');
plot(squeeze(rL2(2,2,:)), squeeze(-rL2(1,2,:)),'m');
plot(squeeze(rL2(2,n,:)), squeeze(-rL2(1,n,:)),'g');
%
plot(squeeze(rL3(2,1,:)), squeeze(-rL3(1,1,:)),'r');
plot(squeeze(rL3(2,2,:)), squeeze(-rL3(1,2,:)),'m');
plot(squeeze(rL3(2,n,:)), squeeze(-rL3(1,n,:)),'g');
%

