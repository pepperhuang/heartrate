function prc_analyze(to_dir)

% Load PRC data
load([to_dir 'phase_data/all_prc.mat'], 'prc_freq', 'prc_scramble');

% PRC is calculated over 18 hours in each direction, 2(18) + 1 = 37
disp_prc = zeros(37, 1);
disp_prc_scramble = zeros(37, 1);

% Putting the PRC data in the appropriate hour bins
disp_prc(8:31) = -prc_freq;
disp_prc(1:7) = -prc_freq(18:24);
disp_prc(32:37) = -prc_freq(1:6);
disp_prc_scramble(8:31) = -prc_scramble;
disp_prc_scramble(1:7) = -prc_scramble(18:24);
disp_prc_scramble(32:37) = -prc_scramble(1:6);

% Plot results
f = figure(1);
set(0, 'CurrentFigure', f);

hold on
p1 = plot((-18):18, disp_prc, 'LineWidth', 4, 'Color', 'blue');
p2 = plot((-18):18, disp_prc_scramble, 'LineWidth', 4, 'Color', 'red');
hold off

% Labels
xlim([-18 18])
xlabel('Relative Time (hours)', 'FontSize', 18)
ylabel('Shift Per Step (minutes)', 'FontSize', 18)
title('Estimated PRC (1 hour window, overall average)', 'FontSize', 24)
legend([p1 p2],{'Real Data', 'Scrambled Activity'}, 'FontSize', 18, 'Location', 'southeast')

% Export results to PDF
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', [to_dir 'phase_data/prc.pdf']);