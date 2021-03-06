% degradation_model.m
% 
% Evaluates probability distributions of cosmogenic exposure dates on a
% degrading moraine.  Based on a conceptual model first described by
% Hallet and Putkonen (1994) and Putkonen and Swanson (2003), but the
% code presented here is entirely original.  Uses an analytical solution 
% for moraine degradation developed by Dr. Nathan Urban, Penn State.  Code 
% written by Patrick Applegate, Penn State (papplegate@psu.edu) and
% modified by Dylan Ward, U Cincinnati (dylan.ward@uc.edu) for episodic
% degradation.
% 
% Code written to run under MATLAB R2008a on an Intel-based Macintosh
% MacBook.  
% 
% Calls the functions m_diffusion.m and ep_diffusion.m, which are provided
% in a separate file.
% 
% This code was written carefully and has been checked for obvious errors.
% However, no warranty of any kind is implied.  The code may not even run
% on your system.  The output from the code should not be trusted without
% testing.  
% 
% Please give proper credit if using this code in research and teaching.
% Derivative works based on this code should include a reference to the
% original paper: 

% Applegate, P. J., Urban, N. M., Laabs, B. J. C., Keller, K., and Alley,
%   R. B.: Modeling the statistical distributions of cosmogenic exposure
%   dates from moraines, Geosci. Model Dev., 3, 293-307,
%   https://doi.org/10.5194/gmd-3-293-2010, 2010.

% Clear all variables, commands, and figures.  Set figures to dock
% automatically.  
clear 
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

% Define parameters that will be tuned during model inversion.  
moraine_age = 20.0;         % ka (10^ 3 yr); true age of moraine
initial_height = 10.0;      % m; initial height of moraine
initial_slope = 25;         % degrees; initial moraine slope angle (Hallet 
                            % and Putkonen (1994) assume 31 degrees; 
                            % Putkonen and Swanson (2003) use 34 degrees; 
                            % 25 degrees may be more reasonable, given 
                            % Putkonen and O'Neal (2006)) 
k = 10^ -3;                 % sq. m/ yr; initial/default topographic diffusion coefficient

% Define variable diffusivity in intervals
% k hist: t_chg(yr) rate(m/yr). Uses default until first defined
% time point.
DO_VARIABLE_DIFF = true;

k_t_mtx = [ 3000,  5e-3;
            4000,  1e-3;
            6000,  5e-3;
            8000,  1e-3 ];
        

% Define other geomorphic parameters that are not part of the inversion.  
erosion_rate = 0;           % mm/ ka; erosion rate of exposed boulders 
boulder_height = .5;       % m; observed height of boulders when sampled
rho_rock = 2.6;             % g/ cm^ 3; density of boulders
rho_till = 2.0;             % g/ cm^ 3; density of till matrix

% Define nuclide production parameters.  
P_spall = 4.97;             % atoms/ g/ yr; surface production rate due to 
                            % spallation (get this from the CRONUS online
                            % calculator described in Balco et al., 2008; 
                            % assumed to be constant over the lifetime of 
                            % the moraine)
P_mu = 0.133;               % atoms/ g/ yr; surface production rate due to 
                            % muons (also get this from the CRONUS
                            % calculator)
decay_const = 4.62* 10^ -7; % yr^ -1; nuclear decay constant of nuclide of
                            % interest (4.62* 10^ -7 for 10Be, following 
                            % Balco et al., 2008, and refs therein)
P_slhl = [5 0.09 ...        % atoms/ g/ yr; sea level, high-latitude 
    0.02 0.02];             % production rates of different cosmic ray 
                            % flux components, following Granger and
                            % Muzikar (2001); for 10Be, about 
                            % [5 0.09 0.02 0.02]
att_length = [160 738 ...   % g/ sq. cm; effective attenuation lengths of 
    2688 4360];             % exponential components of Granger and Muzikar
                            % (2001) production-as-a-function-of-depth
                            % parameterization; for 10Be and 26Al, about
                            % [160 738 2688 4360]

% Define model parameters.  
num_boulders = 10* 10^ 3;    % number of randomly generated synthetic 
                            % boulders (at least 10^ 4; bigger numbers
                            % yield more consistent results, but the model
                            % will take more time to run) 
time_step = 25;             % yr; time step during post-depositional
                            % period (25 yr works well; small values
                            % increase the accuracy of the calculation, 
                            % but also cause the code to run more
                            % slowly)

% Turn plotting on and off.  
plots = 1;                  % If 1, plots the moraine profile, height of 
                            % the moraine's crest as a function of time, 
                            % and a histogram of the modeled exposure 
                            % dates.  If 0, none of these plots are
                            % produced.  
bin_width = .25;              % ka; width of bins in naive age histogram



%% ----- END of INPUTS ----- %%



% Convert all quantities to consistent units.  All lengths should be in
% meters, slopes should be dimensionless, times should be in years, and
% masses should be in grams.  
moraine_age = moraine_age* 10^ 3;               % yr
initial_slope = tand(initial_slope);            % d'less
erosion_rate = erosion_rate* 10^ -6;            % m/ yr
rho_till = rho_till* 100^ 3;                    % g/ cu. m
rho_rock = rho_rock* 100^ 3;                    % g/ cu. m
att_length = att_length* 100^ 2;                % g/ sq. m

% Scale production rates to site.  
P_surf(1) = P_spall; % atoms/ g/ yr
P_surf(2: 4) = P_slhl(2: 4)* P_mu/ sum(P_slhl(2: 4)); 

% Determine length scales for nuclide production.  
L_till = att_length/ rho_till; % m 
L_rock = att_length/ rho_rock; % m

% Set up initial moraine profile
length_step = 1;            % m; space step (1-2 m is best)
length_factor = 3;        % profile length factor (at least 1.5; not 
                            % important, except for large, old moraines) 

% Establish plotting variables times and distances.  
L = initial_height/ initial_slope; % m; half-width of the moraine's base
distances = 0: length_step: (length_factor* L); % m

% Calculate initial moraine profile as a function of distance from the
% moraine's crest.
initial_profile = zeros(1, numel(distances)); % m
for count1 = 1: 1: numel(distances)
    initial_profile(count1) = initial_height- (count1- 1)* ...
        length_step* initial_slope;
    if initial_profile(count1) < 0
        initial_profile(count1) = 0;
    end
end

%% time control for varying degradation rate - modulate diffusivity between periods
times = 0:time_step:moraine_age;        % total time sequence
initial_profile_interval = initial_profile;

                          
% for each interval, create a times_interval array that starts at zero and
% goes to the length of the interval

% Determine height of moraine crest as a function of time
% and final moraine profile.  
% run through time interval times

if DO_VARIABLE_DIFF == true
    num_intervals = size(k_t_mtx,1);
else
    num_intervals = 0;
end

for interval = 0:num_intervals    % interval 0 is that before the first matrix point
    if DO_VARIABLE_DIFF == true 
        if interval == 0
            interval_start = 0;
            k_interval = k;
        else
            interval_start = k_t_mtx(interval,1);
            k_interval = k_t_mtx(interval,2);
        end

        if interval < size(k_t_mtx,1)
            interval_end = k_t_mtx(interval+1,1)-time_step;
        else
            interval_end = moraine_age;
        end
    else     
        interval_start = 0;
        interval_end = moraine_age;
        k_interval = k;
    end

    if interval_end > moraine_age; interval_end = moraine_age; end
    if interval_start > moraine_age; break; end  % stops if we get to the moraine age before the matrix runs out of intervals
    
    interval_length = interval_end - interval_start;
    
    [times_interval, crest_height_interval, final_profile_interval] = ep_diffusion(distances, initial_profile_interval, k_interval, ...
        interval_length, time_step);

    % fill out the crest height over the correct times in the main time
    % matrix
    interval_st_idx = floor(interval_start/time_step) + 1;
    interval_end_idx = interval_st_idx+length(crest_height_interval)-1;
    crest_height(interval_st_idx:interval_end_idx) = crest_height_interval;
    
    % new profile is the start for next interval
    initial_profile_interval = final_profile_interval;

end % for interval

% after time iterations
final_profile = final_profile_interval;

% Establish the initial depth for each boulder. 
final_height = min(crest_height); % m 
max_depth = initial_height- final_height- boulder_height; % m 
initial_depth = max_depth* rand(1, num_boulders); % m 

% Determine the thickness of the erodible shell on each boulder.  This
% thickness depends on the time that each boulder's upper surface is 
% higher than the crest of the moraine, and on the erosion rate.  
shell_thick = zeros(1, num_boulders); % m
if erosion_rate > 0  % don't do these steps if erosion is nil
    for count1 = 1: 1: num_boulders 
        boulder_top = initial_height- initial_depth(count1); % m 
        yn = 0; 
        count2 = 1; 
        while yn == 0 
            if crest_height(count2) <= boulder_top 
                exposure_time = moraine_age- times(count2); % yr
                shell_thick(count1) = exposure_time* erosion_rate; 
                yn = 1; 
            end
            count2 = count2+ 1; 
        end
    end
    initial_shell_thick = shell_thick; % m
end

% Step through time, tracking the nuclide concentration in each boulder.  
boulder_conc = zeros(1, num_boulders); % atoms/ g
for count1 = 2: 1: numel(times) 
    disp(['Calculating time step #', num2str(count1- 1), ' of ',...
        num2str(numel(times)- 1), '... '])
    % Increment concentrations for nuclear decay.  
    boulder_conc = boulder_conc.* exp(-decay_const* time_step); 
    % Step through the list of boulders.  
    for count2 = 1: 1: num_boulders; 
        depth = crest_height(count1)- ...
            (initial_height- initial_depth(count2)); % m 
        % If the boulder is at the surface, 
        if depth <= 0
            P_sample = P_surf.* exp(-shell_thick(count2)./ L_rock); 
            shell_thick(count2) = shell_thick(count2)- ...
                erosion_rate* time_step; 
        % Otherwise, 
        else 
            P_till = P_surf.* exp(-depth./ L_till); 
            P_sample = P_till.* exp(-shell_thick(count2)./ L_rock); 
        end
        % Increment the concentration in the boulder by the production rate
        % during this time step.  
        boulder_conc(count2) = boulder_conc(count2)+ ...
            sum(P_sample)* time_step; 
    end
end

% Calculate the apparent exposure time for each boulder.  
naive_age = -decay_const^ -1* ...
    log(1- ((boulder_conc.* decay_const)./ (P_spall+ P_mu))); % yr

%% Boulder loss as a function of how far exhumed they became

final_bldr_depth = final_height- (initial_height- initial_depth); % m 
% If the boulder is more than one boulder height above the surface,
% it has an increasing probability of loss from the crest
P_remain = ones(size(final_bldr_depth));
P_remain(final_bldr_depth <= -boulder_height) = exp( (final_bldr_depth(final_bldr_depth <= -boulder_height)+boulder_height)/boulder_height);
% figure(2)
% plot(final_bldr_depth,P_remain,'k.')
lost_boulders = P_remain < rand(1, num_boulders);

%% -------- PLOTS -------- %%

% For ease of plotting, convert variables with a time dimension to ka
% (10^3 yr).  
times = times/ 10^ 3; 
naive_age = naive_age/ 10^ 3; 
moraine_age = moraine_age/ 10^ 3; 

if plots == 1
    % Plot the initial (dotted) and final (solid) moraine profiles.
    figure(1)
    subplot(2,2,1)
    plot(distances, initial_profile, 'k--', 'LineWidth', 1.5)
    axis square
    hold on
    plot(distances, final_profile, 'k', 'LineWidth', 1.5)
    xlabel('Distance from moraine crest (m)', 'FontSize', 16, ...
        'FontWeight', 'bold')
    ylabel('Height (m)', 'FontSize', 16, ...
        'FontWeight', 'bold')
    h_leg = legend('Initial profile', 'Final profile');
    legend('boxoff')
    set(h_leg, 'FontSize', 14)
    set(gca, 'FontSize', 14)
    set(gca, 'LineWidth', 1)
%    set(gca, 'Box', 'off')

    % Plot moraine height as a function of time.
    figure(1)
    subplot(2,2,3)
    plot(fliplr(times), crest_height, 'k', 'LineWidth', 1.5)
    axis square
    xlabel('Moraine time (ka)', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('Crest height (m)', 'FontSize', 16, ...
        'FontWeight', 'bold')
    set(gca, 'FontSize', 14, 'XDir','reverse')
    set(gca, 'LineWidth', 1)
%    set(gca, 'Box', 'off')
    
    % Histogram the apparent ages given by the modeled boulders.  
    figure(1)
    subplot(2,2,2)
    nbins = ceil((max(naive_age)- min(naive_age))/ bin_width); 
    h1 = histogram(naive_age, nbins);
    set(h1, 'FaceColor', 'b', 'EdgeColor', 'k')
    hold on
    h1 = histogram(naive_age(~lost_boulders), nbins);
    set(h1, 'FaceColor', 'r', 'EdgeColor', 'k')
    
    ylimits = get(gca, 'YLim'); 
    axis square

    plot([moraine_age moraine_age], ylimits, 'k--', 'LineWidth', 1.5)
    xlabel('Apparent age (ka)', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('Number of boulders', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize', 14)
    set(gca, 'LineWidth', 1)
    set(gca, 'XTickMode', 'auto')
%    set(gca, 'Box', 'off')

  % CDF of the apparent ages given by the modeled boulders.  
    figure(1)
    subplot(2,2,4)
    ct = 0;
    for a = 0:time_step/1000:moraine_age
        ct = ct+1;
        cumdist(ct) = sum(naive_age <= a);
        cumdist_remain(ct) = sum(naive_age(~lost_boulders) <= a);
    end
    plot(0:time_step/1000:moraine_age,cumdist/num_boulders, 'b', 'LineWidth', 1.5)
    hold on
    plot(0:time_step/1000:moraine_age,cumdist_remain/sum(~lost_boulders), 'r', 'LineWidth', 1.5)
    axis square
    xlim([floor(min(naive_age)) moraine_age])
    xlabel('Apparent age (ka)', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('P <= age', 'FontSize', 16, ...
        'FontWeight', 'bold')
    set(gca, 'FontSize', 14)
    set(gca, 'LineWidth', 1)
%    set(gca, 'Box', 'off')
end

beep
