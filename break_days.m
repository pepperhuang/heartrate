function break_days(my_id, from_dir, to_dir)

% PARAMETERS
% Minutes of consecutive wake before sleep is ended
consec_wake = 120;

% Minutes of consecutive sleep before sleep is started
consec_sleep = 120;

% Gaps in sleep (minutes) which are still considered consecutive
sleep_jitter = 5;

% If this many minutes of HR data are missing, allow days to break here
missing_break_cutoff = 120;


% CODE
% File names for HR, steps, sleep data
temp = dir([from_dir my_id '_heartrate*.csv']);
str_hr = temp.name;
temp = dir([from_dir my_id '_minuteSteps*.csv']);
str_steps = temp.name;
temp = dir([from_dir my_id '_minuteSleep*.csv']);
str_sleep = temp.name;

% Open files for reading
f_hr = fopen([from_dir str_hr]);
f_steps = fopen([from_dir str_steps]);
f_sleep = fopen([from_dir str_sleep]);
fgetl(f_hr);
fgetl(f_steps);
fgetl(f_sleep);

days_pre = [];
phase = 1;

% With the HR data...
hr_data = textscan(f_hr,'%s %f','Delimiter',',');
hr_data = datenum(hr_data{1}, 'mm/dd/yyyy HH:MM:SS PM') * 24 * 60;

days_addin = [];

% Look for any gaps larger than missing_break_cutoff
% Add the start and endpoint of those gaps to days_addin
for i = 1:(length(hr_data)-1)
    if (hr_data(i+1) - hr_data(i)) > missing_break_cutoff
        days_addin = [days_addin; hr_data(i), hr_data(i+1)];
    end
end

% With the sleep data...
sleep_data = textscan(f_sleep,'%s %f %f','Delimiter',',');
sleep_data = datenum(sleep_data{1}, 'mm/dd/yyyy HH:MM:SS PM') * 24 * 60;

day_end = sleep_data(1) - 1;
day_start = hr_data(1) - 1;

% Go through all of the sleep data in order
for i = max(find(sleep_data>day_start,1),2):length(sleep_data)
    % If currently asleep, look for a gap in sleep indicating wakefulness
    % and store the start and end points of the day. Go to wake phase
    if phase == 0
        if sleep_data(i) - sleep_data(i-1) >= consec_wake
            phase = 1;
            day_start = sleep_data(i-1) + 1;
            day_end = sleep_data(i) - 1;
        end
    else
        % If currently awake, start looking for enough consecutive sleep to
        % go back to the sleep phase. If there is more wakefulness shortly
        % after entering sleep, keep updating the day endpoint until enough
        % consecutive sleep is observed
        if sleep_data(i) - sleep_data(i-1) <= sleep_jitter
            if sleep_data(i) - day_end > consec_sleep
                phase = 0;
                if day_end > day_start
                    days_pre = [days_pre; day_start, day_end];
                end
            end
        else
            day_end = sleep_data(i) - 1;
        end
    end
end

% The day after the final period of sleep is the final day
days_pre = [days_pre; sleep_data(end) + 1, hr_data(end) + 1];

clearvars sleep_data;

% Go back through the recorded days and insert breaks due to missing HR
% data in the appropriate locations
days = [];
j = 1;
i = 0;
while i < size(days_pre,1)
    i = i + 1;
    % Once all of the fake days are inserted, copy the rest of the real
    % days over into the final list
    if j > size(days_addin,1)
        days = [days; days_pre(i,:), 1]; % 1 indicates this is real sleep
    else
        % Otherwise, keep adding real days until finding one that overlaps
        % with a gap in HR data
        if days_pre(i,2) < days_addin(j,1)
            days = [days; days_pre(i,:), 1]; % 1 indicates this is real sleep
        else
            % Use the gap in HR data to break the day into pieces
            days = [days; days_pre(i,1), days_addin(j,1), 0]; % 0 indicates this is just a gap in HR data
            j = j + 1;
            if j <= size(days_addin,1)
                while days_pre(i,2) >= days_addin(j,1)
                    days = [days; days_addin(j-1,2), days_addin(j,1), 0]; % 0 indicates this is just a gap in HR data
                    j = j + 1;
                    if j > size(days_addin,1)
                        break;
                    end
                end
            end
            while days_pre(i,2) < days_addin(j-1,2)
                i = i + 1;
                if i > size(days_pre,1)
                    break;
                end
            end
            if i <= size(days_pre,1)
                days = [days; days_addin(j-1,2), days_pre(i,2), 1]; % 1 indicates this is real sleep
            end
        end
    end
end

% The nights are the midpoints between consecutive days
nights = zeros(size(days,1)-1,2);
for i = 1:length(nights)
    nights(i,1) = (days(i,2) + days(i+1,1))/(24*120);
end

% The second row in nights stores whether the nights are periods of sleep
% or just gaps in HR data
nights(:,2) = days(1:(end-1),3);
days = days(:,1:2);

clearvars days_pre
clearvars days_addin

% Once sleep has been used to create days, find the start and endpoints for
% heart rate corresponding to those days
i = 1;
j = 1;
max_i = size(days, 1);
days_hr = [];
for k = 1:length(hr_data)
    if j == 1
        if hr_data(k) > days(i,j) - 0.5
            j = 2;
            hr_start = k;
        end
    else
        if hr_data(k) > days(i,j) + 0.49
            days_hr = [days_hr; hr_start, k - 1];
            j = 1;
            i = i + 1;
        end
    end
    if i > max_i
        break;
    end
end

clearvars hr_data;

% Now do the same for steps
i = 1;
j = 1;
days_steps = [];
steps_data = textscan(f_steps,'%s %f','Delimiter',',');
steps_data = datenum(steps_data{1}, 'mm/dd/yyyy HH:MM:SS PM') * 24 * 60;
for k = 1:length(steps_data)
    if j == 1
        if steps_data(k) > days(i,j) - 0.5
            j = 2;
            steps_start = k;
        end
    else
        if steps_data(k) > days(i,j) + 0.49
            days_steps = [days_steps; steps_start, k - 1];
            j = 1;
            i = i + 1;
        end
    end
    if i > max_i
        break;
    end
end

% Save results to file
if ~exist([to_dir 'phase_data'],'dir')
    mkdir([to_dir 'phase_data']);
end
save([to_dir 'phase_data/' my_id '_days.mat'],'days','days_hr','days_steps','nights');