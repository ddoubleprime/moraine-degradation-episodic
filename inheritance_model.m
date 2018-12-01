% inheritance_model.m
% 
% Evaluates probability distributions of cosmogenic exposure dates for
% boulders that contain inherited nuclides.  Code written by Patrick
% Applegate, Penn State (papplegate@psu.edu).  
% 
% Code written to run under MATLAB R2008a on an Intel-based Macintosh
% MacBook.  
% 
% This code was written carefully and has been checked for obvious errors.
% However, no warranty of any kind is implied.  The code may not even run
% on your system.  The output from the code should not be trusted without
% testing.  
% 
% Please give proper credit if using this code in research and teaching.
% Derivative works based on this code should include a reference to the
% original paper.  

% Clear all variables, commands, and figures.  Set figures to dock
% automatically.  
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

% Define parameters that will be tuned during model inversion.  
moraine_age = 20.0;         % ka (10^ 3 yr); true age of moraine
max_pre_time = 100.0;       % ka; maximum time that any individual boulder 
                            % had to acquire inherited nuclides
max_pre_depth = 2.0;        % m; maximum depth of sample point on any 
                            % boulder during predepositional exposure time
pre_slope = 0;              % degrees; slope of surface from which boulders
                            % are derived
                            
% Define other geomorphic parameters that are not part of the inversion.  
erosion_rate = 0.0;         % mm/ ka; erosion rate of boulders on moraine
rho_rock = 2.6;             % g/ cm^ 3; density of boulders
rho_over = 2.0;             % g/ cm^ 3; density of material overlying 
                            % boulders during predepositional exposure time

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
P_slhl = [4.97 0.09 ...     % atoms/ g/ yr; sea level, high-latitude 
    0.02 0.02];             % production rates of different cosmic ray 
                            % flux components, following Granger and
                            % Muzikar (2001); for 10Be, about 
                            % [4.97 0.09 0.02 0.02]
att_length = [160 738 ...   % g/ sq. cm; effective attenuation lengths of 
    2688 4360];             % exponential components of Granger and Muzikar
                            % (2001) production-as-a-function-of-depth
                            % parameterization; for 10Be and 26Al, about
                            % [160 738 2688 4360]

% Define model parameters.  
num_boulders = 1* 10^ 5;    % number of randomly generated synthetic 
                            % boulders (at least 10^ 3; bigger numbers
                            % yield more consistent results, but the model
                            % will take more time to run) 

% Turn plotting on and off.  
plots = 1;                  % If 1, plots the moraine profile, height of 
                            % the moraine's crest as a function of time, 
                            % and a histogram of the modeled exposure 
                            % dates.  If 0, none of these plots are
                            % produced.  
bin_width = 10;             % ka; width of bins in naive age histogram

% Convert all quantities to consistent units.  All lengths should be in
% meters, times should be in years, and masses should be in grams.  Note
% that pre_slope should remain in degrees -- do not reduce this angle to
% its slope equivalent.  
moraine_age = moraine_age* 10^ 3;               % yr
max_pre_time = max_pre_time* 10^ 3;             % yr
erosion_rate = erosion_rate* 10^ -6;            % m/ yr
rho_over = rho_over* 100^ 3;                    % g/ cu. m
rho_rock = rho_rock* 100^ 3;                    % g/ cu. m
att_length = att_length* 100^ 2;                % g/ sq. m

% Scale production rates to site.  
P_surf(1) = P_spall; % atoms/ g/ yr
P_surf(2: 4) = P_slhl(2: 4)* P_mu/ sum(P_slhl(2: 4)); 

% Determine length scales for nuclide production.  
L_over = att_length/ rho_over; % m 
L_rock = att_length/ rho_rock; % m

% Determine predepositional exposure time for each boulder, and the depth
% of the sample point on each boulder during that time.  
pre_time = max_pre_time* rand(1, num_boulders); % yr
pre_depth = max_pre_depth* rand(1, num_boulders); % m

% Calculate the final nuclide concentration in each boulder.  The
% production rate parameterization here follows Dunne et al. (1999), using
% the four-component production rate scheme described by Granger and
% Muzikar (2001).  
boulder_conc = zeros(1, num_boulders); 
for count1 = 1: 1: num_boulders; 
    P_pre = sum(P_surf.* (1- 3.6* 10^ -6* pre_slope^ 2.64).* ...
        exp((-pre_depth(count1)./ L_over).* (1+ pre_slope^ 2/ 5000))); 
    C_pre = (P_pre/ decay_const)* ...
        (1- exp(-decay_const* pre_time(count1)));
    boulder_conc(count1) = C_pre* exp(-decay_const* moraine_age)+ ...
        (P_spall+ P_mu)/ (decay_const+ erosion_rate/ L_over(1))* ...
        (1- exp(-moraine_age* (decay_const+ erosion_rate/ L_over(1)))); 
end

% Calculate the apparent exposure time for each boulder.  
naive_age = -decay_const^ -1* ...
    log(1- ((boulder_conc.* decay_const)./ (P_spall+ P_mu))); % yr

% For ease of plotting, convert variables with a time dimension to ka
% (10^3 yr).  
naive_age = naive_age/ 10^ 3; 
moraine_age = moraine_age/ 10^ 3; 
pre_time = pre_time/ 10^ 3; 

if plots == 1; 
    % Histogram the apparent ages given by the modeled boulders.  
    figure
%     bins = floor(min(naive_age)): bin_width: ceil(max(naive_age)); 
%     bar(bins, histc(naive_age, bins), 'histc')
    nbins = ceil((max(naive_age)- min(naive_age))/ bin_width); 
    hist(naive_age, nbins)
    axis square
    xlabel('Apparent age (ka)', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('Number of boulders', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize', 14)
    set(gca, 'LineWidth', 1)
    set(gca, 'XTickMode', 'auto')
%    set(gca, 'Box', 'off')
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', 'b', 'EdgeColor', 'k')
    hold on
    ylimits = get(gca, 'YLim'); 
    plot([moraine_age moraine_age], ylimits, 'k--', 'LineWidth', 1.5)
end

beep