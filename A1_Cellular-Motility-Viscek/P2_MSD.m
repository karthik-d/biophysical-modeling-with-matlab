%{
Calculate the mean squared displacement (MSD) averaged over all
particles in the system for different values of the noise; use η = 0.4, 0.5, 0.6, 0.7 and
0.8, and at fixed packing fraction φ = 0.5. Use a simulation with at least NT = 5 × 104
steps, and neglect data from the first 2500 steps. Make a plot of MSD vs. the time
window size ∆t for these different values for the noise. As in assignment 0, calculate
the diffusion coefficient D and power law α from the plotted data; make a plot of
D and α as a function of noise η.
%}

eta_values = (0.4:0.1:0.8)';

%% Set the values of other parameters.

% Number of particles
N = 100;

% Define Simulation Parameters
r0      = 2.0;              % Zone of interaction, in units of particle diameters
v0      = 0.05;             % Speed of cells
Nsteps  = 50000;             % Number of steps for the simulation
dt      = 0.005;            % Simulation timestep
beta    = 10000;            % Value of beta for the simulation

% Number of plot steps to skip between consle output
Nplot   = Nsteps/25;

% Packing fraction 
phi = 0.5;

for i = 1:size(eta_values, 1)
    eta = eta_values(i);

    % run simulation.
    [xtotal, ytotal, pols] = vicsek(N, phi, r0, v0, dt, eta, beta, Nsteps, Nplot);

    % compute and store MSDs for each particle.
    MSDtotal = zeros(N, Nsteps);
    for j = 1:N
        MSDtotal(j, :) = compute_msd(xtotal(j, :)', ytotal(j, :)');
    end

    % average the MSD across all particle.
    MSD_average = mean(MSDtotal, 1);
end



% helper function to compute the MSD for a given set of X and Y positions for several time points.
function MSD = compute_msd(xs, ys)
    % number of time points
    NT = size(xs, 1);
    MSD = zeros(NT-1, 1);

    % loop over all possible time window sizes.
    for tw = 1:(NT-1)
	
        %% -- Each loop deals with delT = 1, 2, 3, ... and computes one row of the MSD in the end.
        
        % calculate x and y displacements, separated by tw (time window).
        dxs = xs(1+tw:end) - xs(1:end-tw);
        dys = ys(1+tw:end) - ys(1:end-tw);
        
        % take mean over all displacements.
        dispMean = mean( dxs.^2 + dys.^2 );
        
        % store in MSD array
        MSD(tw) = dispMean;
    end
end