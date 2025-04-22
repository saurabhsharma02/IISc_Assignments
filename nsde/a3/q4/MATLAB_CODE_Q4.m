% Poisson equation for the deflection of a bar with a distributed load
% Ac: cross-sectional area
% E: Young's modulus
% u: deflection
% x: distance measured along the bar's length
% P(x): distributed load
% AcEd2u/dx2 = P(x)

% Parameters
A_c = 0.1; % m^2
E = 200e9; % N/m^2
L = 10; % m
P_x = @(x) 100; % N/m
delta_x = 0.5; % m
N = L/delta_x; % number of elements

% Element stiffness matrix and load vector
k_e = A_c*E/delta_x*[1 -1; -1 1];
f_e = [delta_x/2; delta_x/2]*P_x(delta_x/2);

% Global stiffness matrix and load vector
K = zeros(N+1);
F = zeros(N+1,1);
for i = 1:N
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k_e;
    F(i:i+1) = F(i:i+1) + f_e;
    f_e = [delta_x/2; delta_x/2]*P_x(i*delta_x + delta_x/2);
end
K(1,:) = 0; K(1,1) = 1; F(1) = 0; % boundary condition: u(0) = 0
K(end,:) = 0; K(end,end) = 1; F(end) = 0; % boundary condition: u(L) = 0

% Solve the system of linear equations
u = K\F;

% Plot the deflection profile
x = linspace(0,L,N+1);
plot(x,-u,'LineWidth',2,'color','#0072BD'); % negate u to represent downward deflection
xlabel('Distance along the bar (m)','FontSize',14,'FontWeight','bold');
ylabel('Deflection (m)','FontSize',14,'FontWeight','bold');
title('Deflection of a Bar with a Distributed Load','FontSize',18,'FontWeight','bold','Color','#0072BD');
grid on;
ax = gca;
ax.FontSize = 12;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontWeight = 'bold';
ax.Title.FontWeight = 'bold';
ax.GridColor = [0.5 0.5 0.5];
ax.GridLineStyle = '-';
ax.Box = 'on';
ax.LineWidth = 2;
