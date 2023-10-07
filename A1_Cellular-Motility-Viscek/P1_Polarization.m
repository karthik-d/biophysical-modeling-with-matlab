%{
In a new script, calculate the time-averaged polarization Φ as a function of noise
strength η and packing fraction φ; use φ = 0.25,0.375, and 0.5. For the noise strength,
use values between η = 0.5 and η = 0.7, with points separated by 0.01. Make a plot of
Φ vs η for the three different values of φ from simulations run for NT = 10000 time steps.
For your calculation of Φ, only use velocities from the second half of the simulation
to minimize transient effects that occur at thebeginning of the simulations.
%}

% clear workspace
clear;

% close any open figures
close all;

% clear text from command window
clc;

% Set ranges for noise and packing fraction, respectively.
eta_values = (0.5:0.01:0.7)';
phi_values = [0.25; 0.375; 0.5];


%% Set the values of other parameters.

% Number of particles
N = 100;

% Define Simulation Parameters
r0      = 2.0;              % Zone of interaction, in units of particle diameters
v0      = 0.05;             % Speed of cells
Nsteps  = 10000;             % Number of steps for the simulation
dt      = 0.005;            % Simulation timestep
beta    = 10000;            % Value of beta for the simulation

% Number of plot steps to skip between consle output
Nplot   = Nsteps/25;

plot_colors = ['r', 'g', 'b'];

figure(1), clf, hold on, box on;
for i=1:size(phi_values, 1)
    phi = phi_values(i);

	% collect averaged polarization values.
	avg_pols = zeros(size(eta_values, 1));

	for j=1:size(eta_values, 1)
		eta = eta_values(j);

		% run simulation
		[xtotal, ytotal, pols] = vicsek(N,phi,r0,v0,dt,eta,beta,Nsteps,Nplot);

		% collect polarization value.
		avg_pols(j) = mean(pols(ceil(Nsteps/2):end));
    end

	% plot pol vs. eta, for current packing fraction.
    hold on;
	plot(eta_values, avg_pols, "Color", plot_colors(i));
end

hold off;
xlabel('noise strength $(\eta)$','Interpreter','latex');
ylabel('time-averaged polarization $(<\Phi(t)>)$','Interpreter','latex');
legend("pf=0.250", "pf=0.375", "pf=0.500", ".");
ax = gca;
ax.FontSize = 18;

