% Final Project
% find equations of motion for triple pendulum 3 different ways
% SYMBOLIC SUB VERSION
function [tout,zout] = FPwrite()
%%%
n=3; p.n=n; % set number of links
%%
% set constant values (m's, l's, etc. to EOM)
p.m = 3*ones(n,1) + [zeros(n-1,1); 10]; p.l = ones(n,1); p.d = .5*ones(n,1); p.Ig = ones(n,1); p.g=10;
% set common tspan n ICs n options stuff
tspan = linspace(0,6,2000);
%t0 = (pi/2)*ones(n,1)+[0; rand(n-3,1); 2*rand(2,1)]; 
vals = linspace(0,sqrt(2*pi),n)';
t0 = vals.^2 + (pi/2)*ones(n,1);
td0 = zeros(n,1);
z0 = [t0; td0];
options = odeset('relTol',1e-8,'AbsTol',1e-8);
%%
% get the equations of motion
[A,b] = newton_pendulum_derive(n); %returns all symbolic equations
writeODE(A,b,0); % write RHS file

% solve the ODE
[tout,zout] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
%plotStuff(tout,zout,n,p.l);
totalE = energy(tout, zout, p);
%%
% get the equations of motion a second way
[A,b] = lagrangeDerive(n); %returns all symbolic equations
writeODE(A,b,0); % write RHS file

% solve the ODE
[tout2,zout2] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
%plotStuff(tout,zout,n,p.l);
totalE = energy(tout2, zout2, p);
%%
% get the equations of motion a third way
[A,b] = dirtyLagrangeDerive(n); %returns all symbolic equations
writeODE(A,b,1); % write RHS file

z0 = [z0(1:3)+.00001; z0(4:6)]; % let the initial conditions deviate just a smidge!
% solve the ODE
[tout3,zout3] = ode45(@(t,y)nPendODE(t,y,p), tspan, z0, options);
%plotStuff(tout3,zout3,n,p.l);
totalE = energy(tout3, zout3, p);
%% 
% compare
plotTripleStuff(tout,zout,zout2,zout3,n,p.l);


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

figure(); hold on; axis([-n n -n n]);
title('N Pendulum Motion');
% plot trajectories of ends of 1st and second bars
plot(squeeze(rL(2,1,:)), squeeze(-rL(1,1,:)),'b');
plot(squeeze(rL(2,2,:)), squeeze(-rL(1,2,:)),'m');
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

figure(); hold on; axis([-n n -n n]);
title('N Pendulum Motion');
% % plot trajectories of ends of 1st and second bars
% plot(squeeze(rL(2,1,:)), squeeze(-rL(1,1,:)),'b');
% plot(squeeze(rL(2,2,:)), squeeze(-rL(1,2,:)),'m');
% plot(squeeze(rL(2,n,:)), squeeze(-rL(1,n,:)),'g');
% %
% plot(squeeze(rL2(2,1,:)), squeeze(-rL2(1,1,:)),'b');
% plot(squeeze(rL2(2,2,:)), squeeze(-rL2(1,2,:)),'m');
% plot(squeeze(rL2(2,n,:)), squeeze(-rL2(1,n,:)),'g');
%
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
    h = plot(y,-x,'x-g');
    h2 = plot(y2,-x2,'o-b');
    h3 = plot(y3,-x3,'o-y');
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
plot(y,-x,'x-b');
%
x2 = 0;
y2 = 0;
for l = 1:n
    x2 = [x2, rL2(1,l,end)];
    y2 = [y2, rL2(2,l,end)];
end
plot(y2,-x2,'o-b');


