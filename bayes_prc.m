function [prc_freq_raw, prc_freq_scramble_raw, prc_n, prc_n_scramble] = bayes_prc(my_id,from_dir,to_dir,real_bin_size)

% PRC bin size (minutes)
bin_size = 60;

% Parameters for shifting. By default, compare yesterday to today using
% yesterday's phase as a point of reference; shift_l shifts the starting
% points earlier (days), shift_r shifts the ending point later (days),
% anchor_l shifts the reference point earlier (days)
% Leave these at 0 unless trying to look at correlations
shift_l = 0;
shift_r = 0;
anchor_l = 0;

% Steps data is one of two things needed to compute the PRC...
temp = dir([from_dir my_id '_minuteSteps*.csv']);
str_steps = temp.name;
f_steps = fopen([from_dir str_steps]);
fgetl(f_steps);

% Same setup as bayes_hr_stream
load([to_dir 'phase_data/' my_id '_arrange.mat'],'days_arrange');
load([to_dir 'phase_data/' my_id '_days.mat'],'days_steps');

to_max = size(days_steps,1);
mod_max = size(days_arrange,1);
while days_arrange(mod_max,2) > to_max
    mod_max = mod_max - 1;
end

temp = textscan(f_steps,'%s %f','Delimiter',',');
steps_times = datenum(temp{1}, 'mm/dd/yyyy HH:MM:SS PM') * 24 * 60;
steps_values = temp{2};
clearvars temp;
fclose(f_steps);

refdate = (floor(steps_times(1)/(24*60))-1)*24*60;
steps_times = steps_times - refdate;

% Initialize vectors to store PRC
prc_freq_raw = zeros(24*60/bin_size,1);
prc_n = zeros(24*60/bin_size,1);
prc_freq_scramble_raw = zeros(24*60/bin_size,1);
prc_n_scramble = zeros(24*60/bin_size,1);

% The second piece needed to calculate PRC is the phase predictions
load([to_dir, 'phase_data/' my_id '.mat'],'all_fits');

mod_max = min(mod_max, size(all_fits,1));

% Random permutation of the days to use for scrambled PRC
scramble = zeros(mod_max-shift_r,1);
scramble((2+shift_l):(mod_max-shift_r)) = randperm(mod_max-shift_r-shift_l-1) + shift_l + 1;

% Loop over all nights...
for i = (2+shift_l):(mod_max-shift_r)
    disp([num2str(my_id) ' - ' num2str(i-1-shift_l) '/' num2str(mod_max-shift_r-shift_l-1)]);
    
    % Find the estimated phase difference between the two nights
    obs_shift = (all_fits(i+shift_r,3) - all_fits(i-1-shift_l,3))*real_bin_size;
    
    % Ensure shifts are -12 to 12 hours
    while obs_shift > 12*60
        obs_shift = obs_shift - 24*60;
    end
    while obs_shift < -12*60
        obs_shift = obs_shift + 24*60;
    end
    
    % Now load in the steps data corresponding to that night...
    raw_steps_data = zeros(days_steps(days_arrange(i,2),2)-days_steps(days_arrange(i,1),1),2);
    for j = days_arrange(i,1):(days_arrange(i,2)-1)
        raw_steps_data((days_steps(j,1):days_steps(j,2))-days_steps(days_arrange(i,1),1)+1,1) = steps_times(days_steps(j,1):days_steps(j,2));
        raw_steps_data((days_steps(j,1):days_steps(j,2))-days_steps(days_arrange(i,1),1)+1,2) = steps_values(days_steps(j,1):days_steps(j,2));
    end
    raw_steps_data = raw_steps_data(raw_steps_data(:,1)>0.5,:);

    % Same steps processing as in bayes_hr_stream
    left_min = floor(raw_steps_data(1,1)/bin_size);
    right_max = floor(raw_steps_data(end,1)/bin_size);
    period_offset = [left_min right_max-left_min+1]*bin_size;
    step_int = int32(raw_steps_data(:,1));
    steps_new = zeros(period_offset(2),2);
    steps_new(:,1) = period_offset(1) + (0:(period_offset(2)-1));
    steps_new(step_int-period_offset(1)+1,2) = raw_steps_data(:,2);
    clearvars raw_steps_data;
    [step_times,~,idx] = unique(floor(steps_new(:,1)/bin_size),'stable');

    % Average steps measurements in same bin
    val = accumarray(idx,steps_new(:,2),[],@mean);

    % Calculate which of 24 one-hour bins corresponds to each hour of steps
    % in the day's data. Since the PRC uses activity relative to current
    % phase, subtract off ``current phase''
    step_avg = [int32(mod(step_times - round(all_fits(i-1-anchor_l,3)*real_bin_size/bin_size)-1,24*60/bin_size)+1) val];
    
    % Now do the same things for the scrambled order steps
    raw_steps_data = zeros(days_steps(days_arrange(scramble(i),2),2)-days_steps(days_arrange(scramble(i),1),1),2);
    for j = days_arrange(scramble(i),1):(days_arrange(scramble(i),2)-1)
        raw_steps_data((days_steps(j,1):days_steps(j,2))-days_steps(days_arrange(scramble(i),1),1)+1,1) = steps_times(days_steps(j,1):days_steps(j,2));
        raw_steps_data((days_steps(j,1):days_steps(j,2))-days_steps(days_arrange(scramble(i),1),1)+1,2) = steps_values(days_steps(j,1):days_steps(j,2));
    end
    raw_steps_data = raw_steps_data(raw_steps_data(:,1)>0.5,:);
    left_min = floor(raw_steps_data(1,1)/bin_size);
    right_max = floor(raw_steps_data(end,1)/bin_size);
    period_offset = [left_min right_max-left_min+1]*bin_size;
    step_int = int32(raw_steps_data(:,1));
    steps_new = zeros(period_offset(2),2);
    steps_new(:,1) = period_offset(1) + (0:(period_offset(2)-1));
    steps_new(step_int-period_offset(1)+1,2) = raw_steps_data(:,2);
    clearvars raw_steps_data;
    [step_times,~,idx] = unique(floor(steps_new(:,1)/bin_size),'stable');
    val = accumarray(idx,steps_new(:,2),[],@mean);
    step_scramble = [int32(mod(step_times - round(all_fits(i-1-anchor_l,3)*real_bin_size/bin_size)-1,24*60/bin_size)+1) val];
    
    % Now go through the steps bin by bin and add their contribution to the
    % PRC in each hour, also keeping track of sum of squares in each bin
    for j = 1:size(step_avg,1)
        prc_n(step_avg(j,1)) = prc_n(step_avg(j,1)) + (step_avg(j,2)^2);
        prc_freq_raw(step_avg(j,1)) = prc_freq_raw(step_avg(j,1)) + step_avg(j,2)*obs_shift;
    end
    for j = 1:size(step_scramble,1)
        prc_n_scramble(step_scramble(j,1)) = prc_n_scramble(step_scramble(j,1)) + (step_scramble(j,2)^2);
        prc_freq_scramble_raw(step_scramble(j,1)) = prc_freq_scramble_raw(step_scramble(j,1)) + step_scramble(j,2)*obs_shift;
    end
end

% At the end, divide the sums by sum of squares to get slope coefficient
prc_freq = prc_freq_raw ./ prc_n;
prc_scramble = prc_freq_scramble_raw ./ prc_n_scramble;

% Save results to file
save([to_dir 'phase_data/' my_id '_prc.mat'],'prc_freq','prc_scramble');
end