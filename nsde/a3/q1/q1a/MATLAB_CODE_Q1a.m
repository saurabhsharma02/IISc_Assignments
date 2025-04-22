clc; clear all; close;

load U_64.txt
load U_1024.txt

dx1 = 1/(size(U_64,1) - 1);
x1 = 0:dx1:1;
dx2 = 1/(size(U_1024,1) - 1);
x2 = 0:dx2:1;

% plot
figure(1)
plot(x1, U_64, 'LineWidth', 2)
hold on
plot(x2,U_1024, 'LineWidth', 2)
set(gca,'FontSize',18)
legend('N = 64','N = 1024', 'FontSize', 16, 'Location', 'northwest')
xlabel('x', 'FontSize', 20)
ylabel('U', 'FontSize', 20)
title('Comparison of U for N = 64 and N = 1024', 'FontSize', 24)
grid on
grid minor
box on
