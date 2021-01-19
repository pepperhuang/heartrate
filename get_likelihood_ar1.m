function likelihood = get_likelihood_ar1(params, hr_avg, steps, jumps, N, bin_size)
    % Old code for likelihood contribution due to exponential decay of
    % the effect of activity on heart rate
    %for i = 2:N
    %   hr_z(i) = hr_z(i-1) * exp(-params(5)*jumps(i)) + exp(params(5)*((1:(jumps(i)+0.1))-jumps(i)))*steps(offsets1(i-1):offsets(i));
    %end
    
    % Difference between observed data and prediction without noise
    r = hr_avg(:, 2) - (params(1) + params(2) * cos(pi * (hr_avg(:, 1) - params(3)) / (24 * 30 / bin_size)) + params(4) * steps);
    
    % When there are gaps in HR, noise is the convolution of several
    % individual steps of noise. Prepower computes these up to a max gap
    % size of 10 (for computational efficiency). Terms decay fairly quickly
    % so in practice 10 was enough to lose correlation completely
    prepower = ones(N, 1);
    b_max = min(max(jumps(2:end)) - 1, 10);
    for i = 1:b_max
        prepower(jumps > i) = prepower(jumps > i) + params(6) ^ (2 * i);
    end
    prepower(1) = 1;
    
    % Contribution to likelihood due to error at each step
    likelihood = sum((r(2:end) - (params(6) .^ jumps(2:end)) .* r(1:(end - 1))) .^ 2 ./ (prepower(2:end))) / (params(5) * params(5));
    
    % Contribution to likelihood due to relation between consecutive errors
    % (here using AR(1) error model)
    likelihood = likelihood + sum(log(params(5) * params(5) * prepower)) - log(1 - params(6) * params(6)) + r(1) * r(1) * (1 - params(6) * params(6)) / (params(5) * params(5));
end