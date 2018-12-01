function [times, crest_height, distances, initial_profile, ...
    final_profile] = m_diffusion(initial_height, initial_slope, k, ...
    moraine_age, time_step)

% m_diffusion.m
% 
% Calculates the height of a moraine's crest as a function of time, plus
% the final topographic profile of the moraine.  Assumes a "sawtooth"
% initial profile.  Based on an analytical solution developed by Dr. Nathan
% Urban, Penn State.  
% 
% Syntax: [times, crest_height, distances, initial_profile, ...
%    final_profile] = m_diffusion(initial_height, initial_slope, k, ...
%    moraine_age, time_step)
% times, vector of elapsed time values (yr)
% crest_height, height of moraine crest as a function of the values in the
%   vector times (m)
% distances, vector of distance from the moraine crest (m)
% initial_profile, height of moraine as a function of the values in the
%   vector distances (m)
% final_profile, height of moraine as a function of the values in the
%   vector distances (m)
% initial_height, initial height of the moraine (m)
% initial_slope, initial slope of the moraine sides (d'less)
% k, topographic diffusion coefficient (sq. m/ yr)
% moraine_age, assumed age of moraine (yr)
% time_step, interval between calculations of moraine height (yr)

% Set model variables.  
length_step = 2;            % m; space step (1-2 m is best)
length_factor = 1.5;        % profile length factor (at least 1.5; not 
                            % important, except for large, old moraines) 

% Establish plotting variables times and distances.  
L = initial_height/ initial_slope; % m; half-width of the moraine's base
times = 0: time_step: moraine_age; % yr
distances = 0: length_step: (length_factor* L); % m

% Calculate height of moraine as a function of time.  
crest_height = zeros(1, numel(times)); % m
h0 = initial_height; % m  
for count1 = 1: 1: numel(times); 
    t = times(count1); 
    crest_height(count1) = (h0/ L)* ((2* sqrt(k* t)/ sqrt(pi))* ...
        (exp(-L^ 2/ (4* k* t))- 1)+ ...
        L* erf(L/ (2* sqrt(k* t)))); 
end

% Calculate initial moraine profile as a function of distance from the
% moraine's crest.
initial_profile = zeros(1, numel(distances)); % m
for count1 = 1: 1: numel(distances);
    initial_profile(count1) = initial_height- (count1- 1)* ...
        length_step* initial_slope;
    if initial_profile(count1) < 0;
        initial_profile(count1) = 0;
    end
end

% Calculate final moraine profile as a function of distance from the
% moraine's crest.
final_profile = zeros(1, numel(distances)); % m
h0 = initial_height; % m
t = moraine_age; % yr
for count1 = 1: 1: numel(distances);
    x = distances(count1); % m
    z1 = exp(-(L+ x)^ 2/ (4* k* t))- ...
        2* exp(-x^ 2/ (4* k* t))+ ...
        exp(-(L- x)^ 2/ (4* k* t));
    z2 = (L+ x)* erf((L+ x)/ (2* sqrt(k* t)))- ...
        2* x* erf(x/ (2* sqrt(k* t)))+ ...
        (L- x)* erf((L- x)/ (2* sqrt(k* t)));
    final_profile(count1) = (h0/ (2* L))* ...
        ((2* sqrt(k* t)/ sqrt(pi))* z1+ z2);
end