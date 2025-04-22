M = readmatrix("EulerQ4.csv");
plot (M(:,1),M(:,2));

xlabel('x','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Function plot for step size= 0.1")
grid on
grid minor