load Source_value_100.csv
plot (Source_value_100(:,1),Source_value_100(:,2))

hold on
load Source_value_500.csv
plot (Source_value_500(:,1),Source_value_500(:,2))


hold on
load Source_value_1000.csv
plot (Source_value_1000(:,1),Source_value_1000(:,2))


hold on
load Source_value_1500.csv
plot (Source_value_1500(:,1),Source_value_1500(:,2))



grid on
grid minor
legend('S=100', 'S=500', 'S=1000', 'S=1500')
title("Temperature Distribution at different source values")
xlabel('Radius r')
ylabel('Temperature value T')



figure
load Source_value_396.csv
plot (Source_value_396(:,1),Source_value_396(:,2))
grid on
grid minor
title("Temperature Distribution for source value of 396")