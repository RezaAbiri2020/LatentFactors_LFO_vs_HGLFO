function [trial_start_time,trial_end_time] = get_trial_times(anin,Fs)
% [trial_start_time,trial_end_time] = get_trial_times(anin,Fs)
% 
% function looks for pulses in the analog input signal
%   single pulse - trial start
%   double pulse - trial end
% 
% Inputs:
%   anin - analog input channel
%   Fs - sampling rate 
% 
% Outputs:
%   trial_start_time - vector of trial start times (single pulse times)
%   trial_end_time - vector of trial end times (double pulse times)

N = length(anin);
time = (0:N-1)/Fs;

anin = anin > 2.5;
pulseidx = find(diff(anin)<-.5)+1;
pulsetime = time(pulseidx);

i = 1;
ct = 1;
while i<length(pulseidx),
    % within group of pulses
    num_pulses = sum(pulsetime>=pulsetime(i) & pulsetime<=pulsetime(i)+1.5);
    pulse_groups(ct) = num_pulses;
    group_times(ct) = pulsetime(i);
    group_idx(ct) = pulseidx(i);
    ct = ct + 1;
    i = i + num_pulses;
end

trial_start_idx = (pulse_groups==1);
trial_start_time = group_times(trial_start_idx);

trial_end_idx = (pulse_groups==2);
trial_end_time = group_times(trial_end_idx);
