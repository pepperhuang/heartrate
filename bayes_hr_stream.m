function bayes_hr_stream(my_id, from_dir, to_dir, bin_size)

% Get strings for heart rate, steps data
temp = dir([from_dir my_id '_heartrate*.csv']);
str_hr = temp.name;
temp = dir([from_dir my_id '_minuteSteps*.csv']);
str_steps = temp.name;

% Customizeable parameters for TMCMC
% Number of samples per estimate, burn-in ratio (usual MCMC burn-in)
num_samples = 100000;
burnin_ratio = 0.5;

% Credible interval should contain this fraction of samples
confidence = 0.8;

% Read HR data, steps data
f_hr = fopen([from_dir str_hr]);
f_steps = fopen([from_dir str_steps]);
fgetl(f_hr);
fgetl(f_steps);

% The cutoffs of days are loaded from break_days and arrange_days
load([to_dir 'phase_data/' my_id '_arrange.mat'], 'days_arrange', 'nights_trimmed');
load([to_dir 'phase_data/' my_id '_days.mat'], 'days_hr', 'days_steps');

% HR and steps data do not always align. to_max is the smaller of the two
to_max = min(size(days_hr, 1), size(days_steps, 1));

% If there are more days (determined from sleep) than data, only go up to
% the end of the HR or steps data
mod_max = size(days_arrange, 1);
while days_arrange(mod_max, 2) > to_max
    mod_max = mod_max - 1;
end

% Load HR data to memory
temp = textscan(f_hr, '%s %f', 'Delimiter', ',');
hr_times = datenum(temp{1}, 'mm/dd/yyyy HH:MM:SS PM') * 24 * 60;
hr_values = temp{2};
clearvars temp;

% Load steps data to memory
temp = textscan(f_steps,'%s %f','Delimiter',',');
steps_times = datenum(temp{1}, 'mm/dd/yyyy HH:MM:SS PM') * 24 * 60;
steps_values = temp{2};
clearvars temp;

% Close files for reading
fclose(f_hr);
fclose(f_steps);

% Time 0 is whichever of HR data or steps data begins first
refdate = min([hr_times(1), steps_times(1)]);
refdate = (floor(refdate / (24 * 60)) - 1) * 24 * 60;

% Shift all times by time 0
hr_times = hr_times - refdate;
steps_times = steps_times - refdate;

% Initial parameter guess for TMCMC
% (1): Constant vertical shift in HR
% (2): Amplitude of underlying oscillation
% (3): Horizontal shift of underlying oscillation
% (4): Instantaneous effect of one step on HR
% (5): Measurement noise standard deviation
% (6): Error correlation decay rate
optfit = [80, 10, 940/bin_size, 0.2, 10, 0.5];

% Optimal parameter values and uncertainties are stored in these lists
all_fits = [];
all_stds = [];

% Tracking the most recent fitted phases. These are only used to guess
% starting points for the following day, so that convergence ideally
% happens somewhat faster -- should not affect final results
last_ten = ones(10,1) * optfit(3);
last_ten_var = ones(10,1) * 40000;

% Make sure folder for phase data exists
if ~exist([to_dir 'phase_data'],'dir')
    mkdir([to_dir 'phase_data']);
end
if ~exist([to_dir 'phase_data/' my_id],'dir')
    mkdir([to_dir 'phase_data/' my_id]);
end

% Figure reference for plotting days
f = figure(1);

% Now loop over all nights to analyze...
for i = 1:mod_max
    disp([num2str(my_id) ' - ' num2str(i) '/' num2str(mod_max)]);
    
    % Load in the chunks of HR and steps data corresponding to day i
    raw_hr_data = zeros(days_hr(days_arrange(i, 2), 2) - days_hr(days_arrange(i, 1), 1), 2);
    raw_steps_data = zeros(days_steps(days_arrange(i, 2), 2) - days_steps(days_arrange(i, 1), 1), 2);
    for j = days_arrange(i, 1):days_arrange(i, 2)
        raw_hr_data((days_hr(j, 1):days_hr(j, 2))-days_hr(days_arrange(i, 1),1) + 1, 1) = hr_times(days_hr(j, 1):days_hr(j, 2));
        raw_hr_data((days_hr(j, 1):days_hr(j, 2))-days_hr(days_arrange(i, 1),1) + 1, 2) = hr_values(days_hr(j, 1):days_hr(j, 2));
        raw_steps_data((days_steps(j, 1):days_steps(j, 2)) - days_steps(days_arrange(i, 1), 1) + 1, 1) = steps_times(days_steps(j, 1):days_steps(j, 2));
        raw_steps_data((days_steps(j, 1):days_steps(j, 2)) - days_steps(days_arrange(i, 1), 1) + 1, 2) = steps_values(days_steps(j, 1):days_steps(j, 2));
    end
    
    % Remove any gaps in the data
    raw_hr_data = raw_hr_data(raw_hr_data(:, 1) > 0.5, :);
    raw_steps_data = raw_steps_data(raw_steps_data(:, 1) > 0.5, :);
    
    % Find all unique binned times for HR data
    [hr_avg, ~, idx] = unique(floor(raw_hr_data(:, 1) / bin_size), 'stable');

    % Average HR measurements in same bin
    val = accumarray(idx, raw_hr_data(:, 2), [], @mean);

    % Combine unique times and averages
    hr_avg = [hr_avg, val];

    % Find leftmost and rightmost bin
    left_min = floor(min(raw_hr_data(1, 1), raw_steps_data(1, 1)) / bin_size);
    right_max = floor(max(raw_hr_data(end, 1), raw_steps_data(end, 1)) / bin_size);
    clearvars raw_hr_data;
    
    % Total length of interval containing HR data
    period_offset = [left_min, right_max - left_min + 1] * bin_size;
    
    % Precomputed quantities for efficiency: gaps in consecutive measurements
    jumps = [hr_avg(1, 1) + 0.1; hr_avg(2:end, 1) - hr_avg(1:(end - 1), 1)];

    % " " " ": number of measurements after averaging same times
    N = size(hr_avg, 1);

    % Fill in any gaps in step data with zeros
    step_int = int32(raw_steps_data(:, 1));
    steps_new = zeros(period_offset(2), 2);
    steps_new(:, 1) = period_offset(1) + (0:(period_offset(2) - 1));
    steps_new(step_int - period_offset(1) + 1, 2) = raw_steps_data(:, 2);
    clearvars raw_steps_data;

    % Average steps data into bins (including zeros)
    [step_avg, ~, idx] = unique(floor(steps_new(:, 1) / bin_size) , 'stable');
    val = accumarray(idx, steps_new(:, 2), [], @mean);
    step_avg = [step_avg, val];
    
    % We only need to keep the step bins corresponding to HR data
    step_avg = step_avg(hr_avg(:, 1) - step_avg(1, 1) + 1, :);

    % Start by guessing a repeat of the previous day's fit
    params = optfit;
    
    % Specifically for phase (parameter of interest), refine the guess
    % according to the last 10 phases weighted by their certainty
    params(3) = mean(last_ten ./ last_ten_var) / mean(1 ./ last_ten_var);
    
    % Make sure phase is in valid range
    while params(3) > 2880 / bin_size
        params(3) = params(3) - 1440 / bin_size;
    end
    while params(3) < 1440 / bin_size
        params(3) = params(3) + 1440 / bin_size;
    end
    lb = params(3) - 720 / bin_size;
    rb = params(3) + 720 / bin_size;
    
    % Reasonable minimum values to prevent divergence of the likelihood
    params([1 2 4 5 6]) = max(params([1, 2, 4, 5, 6]), [25, 1, 0.1, 1, 0.1]);

    % On the first day, fit data from scratch
    if i == 1
        fun = @(x) -0.5 * get_likelihood_ar1(x, hr_avg, step_avg(:, 2), jumps, N, bin_size) - ((mod(x(3) - XXX + 720 / bin_size, 1440 / bin_size) - 720 / bin_size) ^ 2) / (2 * (YYY ^ 2));
        
    % On future days, likelihood should incorporate previous day (plus
    % Gaussian noise, sd 1 hr) as a prior
    else
        fun = @(x) -0.5 * get_likelihood_ar1(x, hr_avg, step_avg(:,2), jumps, N, bin_size) - ((mod(x(3) - optfit(3) + 720 / bin_size, 1440 / bin_size) - 720 / bin_size) ^ 2) / (2 * ((stdfit / norminv(0.5 * (1 + confidence))) ^ 2 + (60 / bin_size) ^ 2));
    end

    % TMCMC parameters. Number of chains to run in parallel:
    num_walkers = 500;
    
    % The prior should ensure only physical values of the parameters
    fun_pri = @(x) logical(prod(x([1, 2, 4, 5, 6]) > 0)) && (x(6) < 1);
    
    % Start walkers off as random perturbations around the starting guess
    start_locs = repmat(params', 1, num_walkers) .* (1 + 0.1 * randn(6, num_walkers));
    
    % Ensure that starting locations are in valid range
    start_locs(6, :) = min(start_locs(6, :), 0.999);
    start_locs(3, :) = max(min(start_locs(3, :), rb - 0.1) , lb + 0.1);
    
    % Run GWMCMC to sample from the likelihood
    [samples, ~] = gwmcmc_periodic(start_locs, {fun_pri fun}, num_samples, 'BurnIn', burnin_ratio, 'BinSize', bin_size);
    
    % Best fit is posterior mean of all samples
    optfit= mean(samples, [2 3])';
    
    % To avoid periodic boundary effects, phase was allowed to wander the
    % whole real line. At the end, map back to 1 day interval
    optfit(3) = mod(mean(samples(3, :, end), 'all'), 24 * 60 / bin_size);
    
    % To plot the resulting prediction, use the fitted parameters
    hr_pred = optfit(1) + optfit(2) * cos(pi * (hr_avg(:, 1) - optfit(3)) / (24 * 30 / bin_size)) + optfit(4) * step_avg(:, 2);

    % Finding the uncertainty: start by assuming an uncertainty of 6 hours
    % in either direction, then use binary search (10 iterations) to find interval that
    % contains ``confidence'' fraction of all samples
    stdfit = 6 * 60 / bin_size;
    std_jump = 0.5 * stdfit;
    search_space = mod(samples(3, :, :) - optfit(3), 24 * 60 / bin_size);
    search_space = min(abs(search_space(:)), abs(search_space(:) - 24 * 60 / bin_size));
    for j = 1:10
        if mean(search_space<stdfit, 'all') < confidence
            stdfit = stdfit + std_jump;
        else
            stdfit = stdfit - std_jump;
        end
        std_jump = 0.5 * std_jump;
    end
    
    % Plot the day's data and the resulting predictions
    set(0, 'CurrentFigure', f);
    
    % Time values:
    days = bin_size * hr_avg(:, 1) / 60 - 24 * floor(bin_size * hr_avg(1, 1) / (60 * 24));
    
    % Scatter plot of real data and line plot of predicted data
    yyaxis left;
    scatter(days, hr_avg(:, 2) , 40)
    hold on
    plot(days, hr_pred, 'LineWidth', 2);
    hold off
    
    % Labeling the figure
    the_date = datestr((hr_avg(1, 1) * bin_size + refdate) / (24 * 60), 'yy-mm-dd - hh:MM');
    title(['Heart Rate vs. Model Posterior Mode (' the_date ')'], 'FontSize', 24)
    the_date = erase(the_date,':');
    xlabel('Time (hours)', 'FontSize', 18)
    ylabel('Heart Rate (bpm)', 'FontSize', 18)
    set(gcf, 'Position', [0, 0, 1600, 900]);

    % Determining AM or PM
    rhy_min = floor(mod(optfit(3) * bin_size, 60));
    rhy_hr = mod(floor(optfit(3) * bin_size / 60) , 24);
    rhy_ampm = 'pm';
    min_prez = '';
    if rhy_hr > 11
        rhy_ampm = 'am';
    end
    rhy_hr = mod(rhy_hr - 1, 12) + 1;
    if rhy_min < 10
        min_prez = '0';
    end
    xlim([days(1) - 1, days(end) + 1]);
    yl = ylim;
    ylim([20, yl(2)]);
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.05 * (xl(2) - xl(1)), yl(1) + (yl(2) - yl(1)) * 0.95, [num2str(rhy_hr) ':' min_prez num2str(rhy_min) ' ' rhy_ampm], 'FontSize', 14);
    text(xl(1) + 0.03 * (xl(2) - xl(1)), yl(1) + (yl(2) - yl(1)) * 0.92, ['+/- ' num2str(ceil(stdfit * bin_size)) ' mins'], 'FontSize', 14);
    text(xl(1) + 0.77 * (xl(2) - xl(1)), yl(1) + (yl(2) - yl(1)) * 0.95, 'Sleep Midpoint', 'FontSize', 14);
    text(xl(1) + 0.80 * (xl(2) - xl(1)), yl(1) + (yl(2) - yl(1)) * 0.92, datestr(nights_trimmed(i), 'hh:MM PM'), 'FontSize', 14);
    
    % On the right axis, plot steps data
    lefty = find(steps_new(:, 1) > hr_avg(1, 1) * bin_size, 1);
    yyaxis right;
    plot(steps_new(lefty:end, 1) / 60 - 24 * floor(bin_size * hr_avg(1, 1) / (60 * 24)), steps_new(lefty:end, 2));
    ylim([0, 400])
    ylabel('Steps', 'FontSize', 18)
    
    % Save results to PDF
    set(f, 'PaperOrientation', 'landscape');
    set(f, 'PaperUnits', 'normalized');
    set(f, 'PaperPosition', [0, 0, 1, 1]);
    print(gcf, '-dpdf', [to_dir 'phase_data/' my_id '/' the_date '.pdf']);
    
    % Update list of ``last ten phases'' for initial guess
    last_ten(2:10) = last_ten(1:9);
    last_ten_var(2:10) = last_ten_var(2:10);
    last_ten(1) = optfit(3);
    last_ten_var(1) = stdfit ^ 2;
    
    % Add fits and uncertainties to list
    all_fits = [all_fits; optfit];
    all_stds = [all_stds; stdfit];
end

% Save resulting data
save([to_dir 'phase_data/' my_id '.mat'], 'all_fits', 'all_stds');