function [A,b] = newton_pendulum_derive(n)
% this funciton outputs matrices that can be solved at each timestep for
% theta1dd, theta2dd, and theta3dd
%
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


% could require checking over
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
    a(1:3,i) = a(1:3,i) + (d(i)*tdd(i)*et(:,i) - d(i)*(td(i)^2)*er(:,i));
end

% now that all variables are defined
% do n angular momentum balances
% M = sym(zeros(n,1));
% Hd = sym(zeros(n,1));
% for i = n:-1:1
%     for j = 1:i
%         temp1 = cross(d(i)*er(:,i),(m(i)*g*ei));
%         M(j) = M(j) + temp1(3);
%         temp2 = m(i)*cross(d(i)*er(:,i),a(:,i)) + Ig(i)*tdd(i)*ek;
%         Hd(j) = Hd(j) + temp2(3); % take 3rd component instead of dot product with k
%     end
% end

% n AMB's
for i = 1:n
    % for this angular mom balance
    M=0; % reset moment
    H=0; % reset change in ang mom   
    for j = i:n % the ith AMB has link i:n e.g. 2nd has links 2 and 3
        if i>1 % only subtract vectors if we're not at the origin
            rGmod = rG(:,j) - rL(:,i-1); % construct vector from point of angular mom bal to cg we're dealing with ATM
        else
            rGmod = rG(:,j);
        end
        % add up moment and dh/dt contributions from each bar, i:n
        M = M + cross(rGmod,m(j)*g*ei);
        H = H + cross(rGmod,m(j)*a(:,j)) + Ig(j)*tdd(j)*ek; % accel is wrt fixed frame, but r vector is still wrt point of AMB
    end
    AMB(i) = M(3) == H(3);
end

[A,b] = equationsToMatrix(AMB,tdd);
%A = simplify(A); b=simplify(b);
%tdd = A\b;