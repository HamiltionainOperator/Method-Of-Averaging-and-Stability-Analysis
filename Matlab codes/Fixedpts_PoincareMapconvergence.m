clear; clc; close all;

%% parameters
eps   = 0.001;
gamma = -0.2 ;
alpha = 1;
F     = 1;
T     = 2*pi;

%% Duffing ODE
f = @(t,y) [
    y(2);
   -y(1) + eps*(-2*gamma*y(2) - alpha*y(1)^3 + F*cos(t))
];

%% ---------- FIND POINCARE FIXED POINT ----------
Pmap = @(y) poincare_mismatch(y, f, T);
y_guess = [0.4170; 0.7354];

options = optimoptions('fsolve','Display','iter','TolFun',1e-12);
y_star = fsolve(Pmap, y_guess, options);

disp('Poincare fixed point:')
disp(y_star)
%% -----------------------------------------------

%% grid around fixed point
M = 150;        % number of random points
delta = 1e-3;   % radius of neighborhood

points = y_star' + delta*(2*rand(M,2)-1);


%% iterate Poincare map
iters = 50;

figure; hold on;
for k = 1:size(points,1)

    traj = zeros(iters+1,2);
    traj(1,:) = points(k,:);

    y = points(k,:)';

    for n = 1:iters
        [~, Y] = ode45(f, [0 T], y);
        y = Y(end,:)';
        traj(n+1,:) = y';
    end

    plot(traj(:,1), traj(:,2), '-o', 'LineWidth',1);
end

plot(y_star(1), y_star(2), 'ks', 'MarkerSize',10,'LineWidth',2);

xlabel('x');
ylabel('xdot');
title('Iterated Poincare contraction');
axis equal;
grid on;

%% helper function
function F = poincare_mismatch(y, f, T)
    [~, Y] = ode45(f, [0 T], y);
    yT = Y(end,:)';
    F = yT - y;
end



