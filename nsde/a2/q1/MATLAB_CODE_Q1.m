% load data
T = table2array(readtable('temperature_distribution.txt'));

% define grid parameters
Nx = 128;
Ny = 256;
x_start = 2;
x_end = 3;
y_start = 4;
y_end = 6;

% create x and y arrays
x = linspace(x_start, x_end, Nx);
x = x(2:end-1);
y = linspace(y_start, y_end, Ny);

% create meshgrid
[X, Y] = meshgrid(x, y);

% reshape data
A = reshape(T, Nx-2, Ny);
A = A';

% create figure
figure('position', [100, 100, 800, 500]);

% create heatmap
h = pcolor(X, Y, A);
set(h, 'edgecolor', 'none');

% set colormap
cmap = colormap('parula');
colormap(flipud(cmap));

% add colorbar
cb = colorbar;
set(cb, 'fontsize', 14);
ylabel(cb, 'Temperature (\circC)', 'fontsize', 14);

% set axis labels and title
xlabel('x (m)', 'fontsize', 14);
ylabel('y (m)', 'fontsize', 14);
title('Temperature Distribution', 'fontsize', 16);

% set axis limits
xlim([x_start, x_end]);
ylim([y_start, y_end]);

% set tick labels
set(gca, 'fontsize', 12);

% set color limits
caxis([30, 60]);

% add grid lines
grid on;

% add a box around the plot
box on;



