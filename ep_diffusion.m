function [times, crest_height, final_profile] = ep_diffusion(distances, initial_profile, k, ...
    time_interval, dt_external)

% ep_diffusion.m
% 
% Numerically diffuses a moraine profile as a function of time, based on an input
% initial profile. Diffusivity is input each call so it can change for
% different time periods. Profile is input each period as well so can be
% externally modified eg resteepened by an advective erosion process. 
% 
% Syntax: [times, crest_height, distances, final_profile] = ep_diffusion(distances, initial_profile, k, ...
%    time_interval, dt_external)
% times, vector of elapsed time values in increments of dt_ext (yr)
% crest_height, height of moraine crest as a function of the values in the
%   vector times (m)
% final_profile, height of moraine as a function of the values in the
%   vector distances (m)

% distances, vector of distance from the moraine crest (m) 
    % these should be EQUALLY SPACED as dx used by the diffusion code wil
    % be calculated from the second-first interval in distances
% initial_profile, height of moraine as a function of the values in the
%   vector distances (m)
% k, topographic diffusion coefficient (sq. m/ yr)
% time_interval, period over which to iterate the diffusion (yr)
% dt_external, interval between reports of moraine height (yr)

%% Calculate new moraine profile as a function of distance from the
%% moraine's crest.

% Take input profile and numerically diffuse it for input time period
% calculate dx
dx = distances(2)-distances(1); % does not assume a zero distance for first cell

% internal timestep is independent of dt_ext and is defined based on courant criterion
dt = 1; % replace with courant calculation and then round down to nearest dt_ext multiple

% do the diffusion


%% Calculate height of moraine as a function of time.  

% crest_height = max(profile)