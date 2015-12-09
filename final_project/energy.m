function totalE = energy(time, z, p, str)
%this funciton calculates the energy

% unpack stuff
n=p.n;
t = z(:,1:n); td = z(:,n+1:2*n);
l=p.l; m=p.m; d=p.d; Ig=p.Ig; g=p.g;
% define unit vectors
ei = [1 0 0]'; ej = [0 1 0]'; ek = [0 0 1]';
% unit vectors at every timestep
for j = 1:length(time)
    for i = 1:n
        er(1:3,i,j) = [cos(t(j,i)); sin(t(j,i)); 0];
        et(1:3,i,j) = [-sin(t(j,i)); cos(t(j,i)); 0];
    end
end
% velocity of CGs of bars at every timestep and Kinetic Energy
vG = zeros(3, n, length(time)); % initialize
KE = zeros(length(time),1);
for timeInd = 1:length(time)
    % velocity of CGs
    for j=1:n
        vOprime=0;
        for k=1:j-1 % contributions from bars above jth bar
            vOprime = vOprime + td(timeInd,k)*l(k)*et(1:3,k,timeInd);
        end
        % contribution from jth bar
        vG(:,j,timeInd) = vOprime + td(timeInd,j)*d(j)*et(1:3,j,timeInd);
    end
    
    % kinetic energy at every timestep, add contributions from each bar
    for i = 1:n
        KE(timeInd) = KE(timeInd) + .5*m(i)*norm(vG(:,i,timeInd))^2 + .5*Ig(i)*td(timeInd,i)^2;
    end
end

% potential energy
PE_ith = zeros(n,1);
PE = zeros(length(time),1);
for timeInd = 1:length(time)
    for i=1:n % for each bar
        h=0;
        for j=1:i-1 % for each bar above ith bar
            h = h + l(j)*cos(t(timeInd,j)); % contribution from bars above bar of interest
        end
        h=h+d(i)*cos(t(timeInd,i)); % contribution from bar of interest
        PE_ith(i) = -m(i)*g*h; % for the ith bar
    end
    PE(timeInd) = sum(PE_ith);
end

%PE = PE-PE(1); % potential energy is relative
totalE = KE + PE;

%% plot energy
titleStr = sprintf('Energy as a function of time, %s',str);
figure(); hold on;
title(titleStr);
plot(time,KE);
plot(time,PE);
plot(time, totalE);
legend('KE','PE','total');

figure(); hold on;
title(titleStr);
plot(time, totalE);
legend('total');

