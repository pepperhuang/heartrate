function [avg_diffs, abs_diffs, std_list] = phase_sleep(my_id, to_dir, bin_size)

% Load in list of fits from Bayesian UQ
load([to_dir 'phase_data/' my_id '.mat'],'all_fits', 'all_stds');

% Initialize vectors for results
avg_diffs = zeros(1, 3);
abs_diffs = zeros(1, 3);
avg_diffs(1) = -99;

% Calculating average uncertainty and unit converting to hours
std_list = mean(all_stds) * bin_size / 60;

% Assuming at least one night has been successfully fit...
if size(all_fits, 1) > 1
    
    % Load sleep midpoints
    load([to_dir 'phase_data/' my_id '_arrange.mat'], 'nights_trimmed');
    
    % Some nights at the end may not have been fit due to HR or step data
    % ending earlier than sleep data
    nights_trimmed = nights_trimmed(1:size(all_fits, 1), :);
    
    % To find time of day, subtract the floor (units of days)
    sleep_mids = (nights_trimmed - floor(nights_trimmed));
    
    % Now extract phases of all the fits, unit convert to days
    phases = mod(all_fits(:, 3) / (60 / bin_size) - 12, 24) / 24;
    
    % Start making plot of sleep midpoint vs. phase estimate
    f = figure(1);
    
    % Plotting the sleep midpoints
    scatter(nights_trimmed,sleep_mids, 30, [0.7 0.7 0.7], 'filled', 'LineWidth', 1);
    
    % Plotting the estimated phases, color coded by uncertainty.
    % Uncertainty cutoffs here are 2 hrs and 3 hrs
    hold on;
    scatter(nights_trimmed(all_stds > 3 * 60 / bin_size), phases(all_stds > 3 * 60 / bin_size), 40, [0.9 0.1 0], 'x', 'LineWidth', 1);
    scatter(nights_trimmed(and(all_stds >= 2 * 60 / bin_size, all_stds <= 3 * 60 / bin_size)), phases(and(all_stds >= 2 * 60 / bin_size, all_stds <= 3 * 60 / bin_size)), 60, [0.8 0.8 0], 'x', 'LineWidth', 2);
    scatter(nights_trimmed(all_stds < 2 * 60 / bin_size), phases(all_stds < 2 * 60 / bin_size), 80, [0 0.8 0.3], 'x', 'LineWidth', 3);
    hold off;
    
    % Calculate the average difference between sleep midpoint and predicted
    % phase in a few different ways: (1) all points, (2) only points with
    % uncertainty < 3 hours, (3) only points with uncertainty < 2 hours
    avg_diffs(1) = mean(mod(phases - sleep_mids + 0.5, 1) - 0.5) * 24;
    avg_diffs(2) = mean(mod(phases(all_stds < 3 * 60 / bin_size) - sleep_mids(all_stds < 3 * 60 / bin_size) + 0.5, 1) - 0.5) * 24;
    avg_diffs(3) = mean(mod(phases(all_stds < 2 * 60 / bin_size) - sleep_mids(all_stds < 2 * 60 / bin_size) + 0.5, 1) - 0.5) * 24;
    
    % Same for average absolute difference
    abs_diffs(1) = mean(abs(mod(phases - sleep_mids + 0.5, 1) - 0.5)) * 24;
    abs_diffs(2) = mean(abs(mod(phases(all_stds < 3 * 60 / bin_size) - sleep_mids(all_stds < 3 * 60 / bin_size) + 0.5, 1) - 0.5)) * 24;
    abs_diffs(3) = mean(abs(mod(phases(all_stds < 2 * 60 / bin_size) - sleep_mids(all_stds < 2 * 60 / bin_size) + 0.5, 1) - 0.5)) * 24;
    
    % Setting axes and labeling the plot
    ylim([0 1])
    xlim([nights_trimmed(1) nights_trimmed(end)])
    cur_month = month(datetime(nights_trimmed(1), 'ConvertFrom', 'datenum'));
    cur_year = year(datetime(nights_trimmed(1), 'ConvertFrom', 'datenum'));
    mon_lab = [];
    if day(datetime(nights_trimmed(1), 'ConvertFrom', 'datenum')) == 1
        mon_lab = datetime(cur_year, cur_month, 1);
    end
    while 1
        cur_month = cur_month + 1;
        if cur_month > 12
            cur_month = 1;
            cur_year = cur_year + 1;
        end
        newdate = datetime(cur_year, cur_month, 1);
        if newdate <= datetime(nights_trimmed(end), 'ConvertFrom', 'datenum')
            mon_lab = [mon_lab newdate];
        else
            break;
        end
    end
    set(gca, 'XTick', datenum(mon_lab));
    datetick('x', 'mmm ''yy', 'keeplimits', 'keepticks')
    datetick('y', 'HH pm', 'keeplimits')
    xtickangle(90)
    xlabel('Date', 'FontSize', 18)
    ylabel('Time', 'FontSize', 18)
    title(['Sleep Midpoint vs. Estimated Phase (' num2str(my_id) ')'], 'FontSize', 24)
    
    % Save figure file
    savefig([to_dir 'phase_data/' my_id '.fig']);
    
    % Also export result to PDF
    set(f,'PaperOrientation','landscape');
    set(f,'PaperUnits','normalized');
    set(f,'PaperPosition', [0, 0, 1, 1]);
    print(gcf, '-dpdf', [to_dir 'phase_data/' num2str(my_id) '.pdf']);
end