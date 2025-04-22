%% u(x) Plot at t=0 time
fplot (@(x) sin(4*x)+sin(x));
grid on;
grid minor;
xlim([0 2*pi]);
xlabel('x');
ylabel('u(x)');
title("Plot of u(x,0)");


figure
%% Error plot
load('solution.dat');
rd1 = solution(1,1);
rd2 = solution(5,1);

N = solution(1:4,2);

E1 = solution(1:4,3);
E2 = solution(5:8,3);

loglog(N, E1);
hold on;
loglog(N, E2);
xlabel('N', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
legend('rd=0.5', 'rd=0.166667');
title('Error Plot')

figure
%%
% Load the solution data
load('solution.dat');
rd1 = solution(1,1);
rd2 = solution(5,1);

% Extract the values of N and E for the two different rd
N = solution(1:4,2);
E1 = solution(1:4,3);
E2 = solution(5:8,3);

% Plot the errors as a function of N
loglog(N, E1, 'bo-', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
loglog(N, E2, 'rs--', 'LineWidth', 2, 'MarkerSize', 10);

% Add lines representing the order of accuracy for both radii
slope1 = log(E1(end)/E1(1)) / log(N(end)/N(1));
slope2 = log(E2(end)/E2(1)) / log(N(end)/N(1));
order1 = slope1 * log10(N) + log10(E1(1));
order2 = slope2 * log10(N) + log10(E2(1));
loglog(N, 10.^order1, 'k--', 'LineWidth', 2);
loglog(N, 10.^order2, 'k-.', 'LineWidth', 2);

% Add axis labels and a legend
xlabel('N', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
legend('r_d = 0.5', 'r_d = 0.1667', 'Order of accuracy (r_d = 0.5)', 'Order of accuracy (r_d = 0.1667)', 'FontSize', 14);


