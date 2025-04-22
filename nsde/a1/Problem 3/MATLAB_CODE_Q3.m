load criticallydamped.txt
plot (criticallydamped(:,1), criticallydamped(:,2));

hold on
load overdamped.txt
plot (overdamped(:,1), overdamped(:,2));

hold on
load underdamped.txt
plot (underdamped(:,1), underdamped(:,2));
xlabel('Time t (sec)','FontWeight','bold')
ylabel('Displacement x (m)','FontWeight','bold')
grid on
grid minor
legend('Critically-damped(c=40)','Over-damped(c=200)','Under-damped(c=5)');

figure
load convergence.dat
plot (convergence(:,1), convergence(:,2));