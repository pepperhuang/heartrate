%% PREPROCESSING
% Location of input files
from_dir = 'G:/Clark/Desktop/2018 Input/';

% Location to place output files
to_dir = 'G:/Clark/Desktop/2018 Output/';

% Bin size in minutes
bin_size = 5;

% Find all .csv files in input; put unique filenames in file_ids
all_files = dir([from_dir '*.csv']);
file_ids = zeros(1, length(all_files));
for i = 1:length(all_files)
    j = strfind(all_files(i).name, '_');
    if ~isempty(j)
        j = j(1) - 1;
        file_ids(i) = str2double(all_files(i).name(1:j));
    end
end
file_ids = unique(file_ids(file_ids > 0));
if ~exist(to_dir, 'dir')
    mkdir(to_dir);
end

%% MAIN LOOP

% Vectors for list of (1) average difference between estimated phase and
% sleep midpoint, (2) average absolute difference between estimated phase
% and sleep midpoint, and (3) average uncertainty (confidence interval width) for each subject
avg_diffs = [];
abs_diffs = [];
std_list = [];

% Loop over file IDs...
for i = 1:length(file_ids)
    disp(['Loading ID ' num2str(file_ids(i))]);
    
    % Script 1/3: use sleep data to break steps and heart rate into days
    try
        break_days(num2str(file_ids(i)), from_dir, to_dir);
    catch
        disp("Error: break_days")
    end
    
    % Script 2/3: ``days'' as determined by the previous script and real
    % time days may not line up perfect. Index ``days'' according to the
    % real life day correspondence
    try
        arrange_days(num2str(file_ids(i)), from_dir, to_dir);
    catch
        disp("Error: arrange_days")
    end
    
    % Script 3/3: Bayesian uncertainty quantification
    try
        % Estimate phase for all nights of sleep
        bayes_hr_stream(num2str(file_ids(i)), from_dir, to_dir, bin_size);
        
        % Plots of phase vs. sleep midpoint over time
        [new_diffs, new_diffs_2, new_std] = phase_sleep(num2str(file_ids(i)), to_dir, bin_size);
        
        % If at least one day of data was analyzed...
        if new_diffs(1) > -99
            
            % Add new results to memory
            avg_diffs = [avg_diffs; new_diffs];
            abs_diffs = [abs_diffs; new_diffs_2];
            std_list = [std_list; new_std];
        end
    catch
        disp("Error: bayes_hr_stream")
    end
     
end

%% HISTOGRAM OF FIT RESULTS COMPARED TO SLEEP MIDPOINT

f = figure(3);
set(0, 'CurrentFigure', f);
figuresize(30, 10)
subplot(1, 3, 1)
histogram(avg_diffs(:, 1), -12:0.5:12, 'FaceColor', 'black', 'FaceAlpha', 1)
xlabel('Avg Residual (hours)', 'FontSize', 18)
xlim([-12 12])
pbaspect([1 1 1])
xbd = xlim;
ybd = ylim;
text(xbd(2) - 0.15 * (xbd(2) - xbd(1)), ybd(2) - 0.07 * (ybd(2) - ybd(1)), "(A)", 'FontSize', 16)
subplot(1, 3, 2)
histogram(abs_diffs(:, 1), 0:0.25:12, 'FaceColor', 'black', 'FaceAlpha', 1)
xlabel('Avg Abs Residual (hours)', 'FontSize', 18)
xlim([0 12])
pbaspect([1 1 1])
xbd = xlim;
ybd = ylim;
text(xbd(2) - 0.15 * (xbd(2) - xbd(1)), ybd(2) - 0.07 * (ybd(2) - ybd(1)), "(B)", 'FontSize', 16)
subplot(1, 3, 3)
histogram(std_list, 0:0.25:12, 'FaceColor', 'black', 'FaceAlpha', 1)
pbaspect([1 1 1])
xlabel('Avg Error (hours)', 'FontSize', 18)
xlim([0 12])
xbd = xlim;
ybd = ylim;
text(xbd(2) - 0.15 * (xbd(2) - xbd(1)), ybd(2) - 0.07 * (ybd(2) - ybd(1)), "(C)", 'FontSize', 16)
print(gcf, '-dpdf', 'fig3.pdf');

%% USING DATA TO GENERATE THE GLOBAL PHASE RESPONSE CURVE

first_go = true;

% Loop over all IDs...
for i = 1:length(file_ids)
    disp(['Loading ID ' num2str(file_ids(i))]);
    load([to_dir, 'phase_data/' num2str(file_ids(i)) '.mat'], 'all_fits');
    if size(all_fits, 1) < 50
        disp('-- Subject has less than 50 nights of data.')
    else
        
        % For subjects with at least 50 nights of data, calculate PRC
        [a, b, c, d] = bayes_prc(num2str(file_ids(i)), from_dir, to_dir, bin_size);
        
        % prc_freq and prc_scramble keep track of the PRC across all
        % subjects. prc_scramble uses step data which has been randomly
        % permuted to verify that the PRC does not rely on activity trends
        if first_go
            prc_freq = a;
            prc_scramble = b;
            prc_n = c;
            prc_n_scramble = d;
            first_go = 0;
        else
            prc_freq = prc_freq + a;
            prc_scramble = prc_scramble + b;
            prc_n = prc_n + c;
            prc_n_scramble = prc_n_scramble + d;
        end
        clearvars a
        clearvars b
        clearvars c
        clearvars d
    end
end

% To get the total PRC, divide total sum by sum of squared steps
prc_freq = prc_freq ./ prc_n;
prc_scramble = prc_scramble ./ prc_n;

% Store resulting PRC data
save([to_dir 'phase_data/all_prc.mat'], 'prc_freq', 'prc_scramble');

% Plot PRC
prc_analyze(to_dir);

%% FITTING PHASE RESPONSE CURVE PARAMETERS

for i = 1:length(file_ids)
    disp(['Loading ID ' num2str(file_ids(i))]);
    load([to_dir, 'phase_data/' num2str(file_ids(i)) '.mat'], 'all_fits');
    if size(all_fits, 1) < 50
        disp('-- Subject has less than 50 nights of data.')
    else
        bayes_prc_fit(num2str(file_ids(i)), from_dir, to_dir);
    end
end