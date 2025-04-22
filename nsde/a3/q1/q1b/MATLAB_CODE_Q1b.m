clc; clear all; close;
% define ticks
set(0,'defaultaxesTickLength',[0.03 0.03]);
load Ua_1024.txt
load Ub_1024.txt

dx = 1/(size(Ua_1024,1) - 1);
x = 0:dx:1;

%plot 
figure (1)
plot(x, Ua_1024, 'LineWidth', 2)
set(gca,'FontSize',18)
xlabel('x', 'FontSize', 20)
ylabel('U', 'FontSize', 20)
title('Variation of U for N = 1024 and alpha =0.001', 'FontSize', 24)
grid on
grid minor
box on


% plot
figure(2)
plot(x, Ua_1024, 'LineWidth', 2)
hold on
plot(x,Ub_1024, 'LineWidth', 2)
set(gca,'FontSize',18)
legend('Part (a)','Part (b)', 'FontSize', 16, 'Location', 'northwest')
xlabel('x', 'FontSize', 20)
ylabel('U', 'FontSize', 20)
title('Comparison of U for N = 1024 for alpha = 0 (part a) and alpha =0.001 (part b)', 'FontSize', 24)
grid on
grid minor
box on
