% Assignment-3; Question-1
% Explicit; Parabolic
alpha = [0; 0.001];
N = [64; 1024]; % Number of gridpoints
nStep = [10; 50; 100; 150; 188]; % Plottable data available at these timesteps

str_legend = cell(length(nStep),1); % cell array to store legend entries
for a = 1:length(alpha)
    for i = 1:length(N)
        figure
        for j = 1:length(nStep)
            fileID1 = fopen(sprintf('alpha_%1.3f_u_N_%1.0f_nStep_%1.0f.txt',...
                alpha(a), N(i), nStep(j)), 'r');
            u = fscanf(fileID1, '%f');
            fileID2 = fopen(sprintf('position_N_%1.0f.txt', N(i)), 'r');
            pos = fscanf(fileID2, '%f');

            str_legend{j} = sprintf('nStep = %1.0f',nStep(j));
            plot(pos, u)
            hold on
        end
        grid on
        grid minor
        ylim([-3 3])
        title(sprintf('Velocity vs Position for \\alpha = %1.3f, at N = %1.0f', alpha(a), N(i)))
        xlabel('Position, x[m]')
        ylabel('Velocity, [m/s]')
        legend(str_legend, 'Position', [0.1 -0.28 1 1])
    end
end