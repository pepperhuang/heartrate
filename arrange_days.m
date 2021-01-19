function arrange_days(my_id, ~, to_dir)

% PARAMETERS
% How many minutes (minimum) to use on each side of sleep for estimation
min_separation = 18 * 60;

% CODE
% Load days/nights from break_days.m
load([to_dir 'phase_data/' my_id '_days.mat'], 'days', 'nights');

% Days_arrange will store the start/end day surrounding periods of sleep
days_arrange = zeros(size(days, 1) - 1, 2);

% Initialize as including exactly one day on each side
days_arrange(:, 1) = 1:size(days_arrange, 1);
days_arrange(:, 2) = 1 + days_arrange(:, 1);

% For each period of sleep...
for i = 1:size(days_arrange, 1)
    my_midpt = (days(i+1, 1) + days(i, 2)) / 2;
    
    % Include more days to the left until enough data are present
    while days_arrange(i, 1) > 1
        if days(days_arrange(i, 1) - 1, 2) > my_midpt - min_separation
            days_arrange(i, 1) = days_arrange(i, 1) - 1;
        else
            break;
        end
    end
    
    % Include more days to the right until enough data are present
    while days_arrange(i, 2) < size(days, 1)
        if days(days_arrange(i, 2) + 1, 1) < my_midpt + min_separation
            days_arrange(i, 2) = days_arrange(i, 2) + 1;
        else
            break;
        end
    end
end

% Phase is only estimated centered at sleep (not during breaks in HR);
% now that HR gaps have been used to set start/end indices, can discard
days_arrange = days_arrange(nights(:, 2) > 0.5, :);
nights_trimmed = nights(nights(:, 2) > 0.5, 1);

% If any periods of study strictly include others, discard the smaller
i = 0;
n = size(days_arrange, 1);
while i < n
    i = i + 1;
    if sum(and(days_arrange(:, 1) <= days_arrange(i, 1), days_arrange(:, 2) >= days_arrange(i, 2))) > 1
        days_arrange(i, :) = [];
        nights_trimmed(i) = [];
        n = n - 1;
    end
end

% Store output
save([to_dir 'phase_data/' my_id '_arrange.mat'], 'days_arrange', 'nights_trimmed');