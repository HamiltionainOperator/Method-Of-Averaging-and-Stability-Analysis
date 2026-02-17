clear; clc; close all;

%% parameters
eps = 0.001;
gamma = 0.5;
alpha = 1;
F = 1;
T = 2*pi;

%% Duffing ODE (2D system)
f = @(t,y) [
    y(2);
   -y(1) + eps*(-2*gamma*y(2) - alpha*y(1)^3 + F*cos(t))
];

%% shared time grid (important!)
tspan = linspace(0,1000*T,50000);

%% two nearby initial conditions
y0A = [1;0];
y0B = [1+1e-4;0];   % tiny perturbation

%% integrate both
[tA, YA] = ode45(f, tspan, y0A);
[tB, YB] = ode45(f, tspan, y0B);

%% separation distance
diff = YA - YB;
d = sqrt(sum(diff.^2,2));
logd = log(d);
mask = (tA > 3000) & (tA < 8000);
      
p = polyfit(tA(mask), logd(mask), 1);
lambda = p(1);

fprintf('Estimated decay rate = %.6f\n', lambda);
%% plot separation
figure;
plot(tA, logd);
xlabel('t');
ylabel('log distance');
title('Growth of perturbation');
grid on;
