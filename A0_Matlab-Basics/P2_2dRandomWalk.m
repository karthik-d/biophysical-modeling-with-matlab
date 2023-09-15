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
fileData = textscan(fobj, "%f %f %f");
celldisp(fileData);

% close file object
fclose(fobj);

%{
% parse file data to get variables of interest
t = ;
x = ;
y = ;

% get number of time points in trajectory using length function
NT = length(t);


%% Part 2. Calculate the MSD

% create MSD array (y-axis of MSD plot)
MSD = zeros(NT-1,1);

% loop over the different possible time windows, calculate MSD for each
% time window size
for ii = 1:(NT-1)
	% calculate x displacements, separated by ii indices
	dx = x(1+ii:end) - x(1:end-ii);
	
	% calculate y displacements similarly
	dy = ;
	
	% take mean over all displacements
	dispMean = mean(  );
	
	% store in MSD array
	MSD(ii) = dispMean;
end

% create deltaT array, using a for loop or vectorization
deltaT = ;

%% Part 3. Plot MSD on a log-log sclae

% open figure window
figure(1), clf, hold on, box on;

% plot curve, add units to axes, etc


%% Part 4. Determine alpha and diffusion coefficient

% get coefficients for line of first few MSD points using polyfit



% store slope and diffusion coefficient
alpha = ;
diffusionCoefficient = ;

% print results to console
disp(['alpha is = ' num2str(alpha)]);
disp(['measured diffusion coefficient is = ' num2str(diffusionCoefficient)]);


%% Part 5. Compare expected and measured diffusion coefficient

% store step size and step duration
stepSize = ;
stepDuration = ;

% calculate expected diffusion coefficient
expectedD = stepSize*stepSize/stepDuration;

% print result to console
disp(['expected diffusion coefficient is = ' num2str(expectedD)]);

%}


