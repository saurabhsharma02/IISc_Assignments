%%
% load data
load('solution.dat');
rd1 = solution(1,1);
rd2 = solution(5,1);
N = solution(1:4,2);
E1 = solution(1:4,3);
E2 = solution(5:8,3);

% plot data
loglog(N, E1, 'b-', 'LineWidth', 2);
hold on;
loglog(N, E2, 'r--', 'LineWidth', 2);
legend(sprintf('rd=%.2f', rd1), sprintf('rd=%.6f', rd2), 'Location', 'northeast');
xlabel('log N', 'FontSize', 14);
ylabel('log(E(N))', 'FontSize', 14);
title('N vs E(N) graph', 'FontSize', 16);
grid on;

ax = gca; 
ax.FontSize = 12;
ax.GridAlpha = 0.3; 
ax.GridLineStyle = '-'; 
ax.GridColor = 'k'; 
