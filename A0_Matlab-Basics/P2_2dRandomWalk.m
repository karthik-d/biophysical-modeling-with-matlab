%% Assignment 0, Problem 2 template

%{
    Problem 2. Mean square-displacement of a 2D random walker

    Submit your completed version of this script
    with the format LASTNAME_assignment0_problem1.m
    to canvas by 11:59 PM on the due date.
%}

% clear workspace
clear;

% close any open figures
close all;

% clear text from command window
clc;

%% Part 1. Read in data from file

% create string with name and location of file in your file path (edit as
% necessary)
fileName = 'data/randomWalker.dat';

% load file into file object
fobj = fopen(fileName, 'r');

% read in data using textscan function (watch out if the file has a
% header!)
fileData = textscan(fobj, "%f %f %f", "Delimiter", ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', 1);

% close file object
fclose(fobj);

% parse file data to get variables of interest
t = fileData{1};
x = fileData{2};
y = fileData{3};

% get number of time points in trajectory using length function
NT = length(t);


%% Part 2. Calculate the MSD

% create MSD array (y-axis of MSD plot)
MSD = zeros(NT-1,1);

% loop over the different possible time windows, calculate MSD for each
% time window size
for ii = 1:(NT-1)
	
	%% -- Each loop deals with delT = 1, 2, 3, ... and computes one row of the MSD at the end.
	
	% calculate x displacements, separated by ii indices
	dx = x(1+ii:end) - x(1:end-ii);
	
	% calculate y displacements similarly
	dy = y(1+ii:end) - y(1:end-ii);
	
	% take mean over all displacements
	dispMean = mean( dx.^2 + dy.^2 );
	
	% store in MSD array
	MSD(ii) = dispMean;
end


% create deltaT array, using a for loop or vectorization
deltaT = t(2:end) - t(1)*ones(NT-1, 1);

%% Part 3. Plot MSD on a log-log sclae

% open figure window
figure(1), clf, hold on, box on;

% plot curve, add units to axes, etc
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

xlabel('Lag Time (seconds)');
ylabel('MSD (microns^2)');

plot(deltaT, MSD, 'Color', 'red', 'LineWidth', 1.5);

%% Part 4. Determine alpha and diffusion coefficient

% get coefficients for line of first few MSD points using polyfit
linearCoeffs = polyfit(log10(deltaT), log10(MSD), 1);

% store slope and diffusion coefficient.
alpha = linearCoeffs(1);
diffusionCoefficient = linearCoeffs(2);

% print results to console
disp(['alpha is = ' num2str(alpha)]);
disp(['measured diffusion coefficient is = ' num2str(diffusionCoefficient)]);


%% Part 5. Compare expected and measured diffusion coefficient

% store step size and step duration
stepSize = sqrt((x(2:end)-x(1:(end-1))).^2 + (y(2:end)-y(1:(end-1)).^2));   % euclidean distances to compute all step sizes.
stepDuration = t(2:end) - t(1:(end-1));
disp(stepDuration(1:10));
disp(stepSize(1:10));

% calculate expected diffusion coefficient
expectedD = stepSize.*stepSize./stepDuration;

disp(expectedD(1:10));
% print result to console
disp(['expected diffusion coefficient is = ' num2str(expectedD)]);

%}


