%% Assignment 1 exmaple script

%{
    An example script to show you how to
    run the vicsek function, which contains
    the velocity-update function
    vicsekvelocities.m

    NOTE: you cannot run this script without
    first completing the function 
    vicsekvelocities

    Once vicsekvelocities is completed, you can run this script to 
    plot the polarization as a function of time during a simulation

    Play around with eta and phi to see how these parameters change the
    polarization over time; also notice what is happening at EARLY times...
%}

% clear workspace
clear;

% close any open figures
close all;

% clear text from command window
clc;

%% Set the values of the parameters

% Number of particles
N = 100;

% Define Simulation Parameters
r0      = 2.0;              % Zone of interaction, in units of particle diameters
v0      = 0.05;             % Speed of cells
Nsteps  = 5000;             % Number of steps for the simulation
dt      = 0.005;            % Simulation timestep
beta    = 10000;            % Value of beta for the simulation

% Number of plot steps to skip between consle output
Nplot   = Nsteps/25;

% Noise parameter
eta = 0.3;

% Packing fraction 
phi = 0.25;

% Run simulation
[xtotal, ytotal, pols] = vicsek(N,phi,r0,v0,dt,eta,beta,Nsteps,Nplot);

% Plot polarization as a function of time
figure(1), clf, hold on, box on;
plot(1:Nsteps,pols,'k-','linewidth',1.75);
xlabel('time steps','Interpreter','latex');
ylabel('$\Phi(t)$','Interpreter','latex');
ax = gca;
ax.FontSize = 18;