function thresh_increment = calc_thresh_increment(thresh_leap, thresh_decay, ...
    refractoryperiod, relrefperiod, tstep)
%thresh_increment calculate the threshold increment that occurs when a
%spike happens
% parameters
% thresh_leap: initial leap in threshold
% thresh_decay: rate of decay of increased threshold
% refractoryperiod: period during which no spikes can occur (infinite
% threshold)
% refractoryperiod: period during which threshold is raised.
% tstep: timestep

% LSS 8 8 2024.
incremented_length = floor((refractoryperiod + relrefperiod)/tstep) ;
incremented_infinity = floor(refractoryperiod/tstep);
% initialise
thresh_increment = zeros([1 incremented_length]) ; 
thresh_increment(1:incremented_infinity) = inf;
thresh_increment(incremented_infinity+1:incremented_length) = ...
    thresh_leap * exp(-thresh_decay * ((1:(incremented_length -incremented_infinity)) * tstep));

end