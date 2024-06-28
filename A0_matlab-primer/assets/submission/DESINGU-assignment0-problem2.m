%% Assignment 0, Problem 2 template

%% PLEASE READ ADDITIONAL SECTION %%

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

plot(deltaT, MSD, 'Color', 'red', 'LineWidth', 1);
hold on;

%% Part 4. Determine alpha and diffusion coefficient


%% (ADDITIONAL begins) ------------------------------------------------------------------------------
%% Since the amount of statistical error increases at high values of delta T,
%% the results in this experiment are illustrated for two cases.
%% - One involving all 5000 data points.
%% - The second involving only the first 600 walks.
%% The experimental and theoretical values of diffusion coefficient are closer for the latter case.

N = 600;
MSD2 = zeros(N-1, 1);
deltaT2 = t(2:N) - t(1)*ones(N-1, 1);

% loop over the different possible time windows, calculate MSD for each
% time window size
for ii = 1:(N-1)
	
	%% -- Each loop deals with delT = 1, 2, 3, ... and computes one row of the MSD at the end.
	
	% calculate x displacements, separated by ii indices
	dx = x(1+ii:end) - x(1:end-ii);
	
	% calculate y displacements similarly
	dy = y(1+ii:end) - y(1:end-ii);
	
	% take mean over all displacements
	dispMean = mean( dx.^2 + dy.^2 );
	
	% store in MSD array
	MSD2(ii) = dispMean;
end
plot(deltaT2, MSD2, 'Color', 'green', 'LineWidth', 1.5, 'LineStyle', '--');
hold on;

%% (ADDITIONAL ends) ------------------------------------------------------------------------------


% get coefficients for line of first few MSD points using polyfit
linearCoeffs = polyfit(log10(deltaT), log10(MSD), 1);
% for first 600 points
linearCoeffs2 = polyfit(log10(deltaT2), log10(MSD2), 1);

% store slope and diffusion coefficient.
alpha = linearCoeffs(1);
diffusionCoefficient = 10^linearCoeffs(2);
% for first 600 points
alpha2 = linearCoeffs2(1);
diffusionCoefficient2 = 10^linearCoeffs2(2);

plot(deltaT, diffusionCoefficient*(deltaT.^alpha), 'Color', 'blue', 'LineWidth', 1.5, 'LineStyle', ':');
hold on;

plot(deltaT2, diffusionCoefficient2*(deltaT2.^alpha2), 'Color', 'black', 'LineWidth', 1.5, 'LineStyle', ':');
hold off;

legend('Experimental (all walks)', 'Experimental (first 600 walks)', 'Theoretical (all walks)', 'Theoretical (first 600 walks)');
fontsize(22, 'points');
xlabel('log(delT) - unitless');
ylabel('log(MSD) - unitless');

% print results to console
disp(['alpha is = ' num2str(alpha)]);
disp(['alpha is (based on first 600 walks) = ' num2str(alpha2)]);
disp(['measured diffusion coefficient is = ' num2str(diffusionCoefficient)]);
disp(['measured diffusion coefficient is (based on first 600 walks) = ' num2str(diffusionCoefficient2)]);


%% Part 5. Compare expected and measured diffusion coefficient

% store step size and step duration
stepSize = MSD(1);
stepDuration = t(2)-t(1);

% calculate expected diffusion coefficient
expectedD = stepSize*stepSize/stepDuration;

% print result to console
disp(['expected diffusion coefficient is = ' num2str(expectedD)]);
disp(" ");
disp('Since the amount of statistical error increases at high values of delta T, the results in this experiment are illustrated for two cases. One involving all 5000 data points. The second involving only the first 600 walks. The experimental and theoretical values of diffusion coefficient are closer for the latter case.');

%% PLEASE READ ADDITIONAL SECTION %%


%% Assignment 0, Problem 2 template

%% PLEASE READ ADDITIONAL SECTION %%

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

plot(deltaT, MSD, 'Color', 'red', 'LineWidth', 1);
hold on;

%% Part 4. Determine alpha and diffusion coefficient


%% (ADDITIONAL begins) ------------------------------------------------------------------------------
%% Since the amount of statistical error increases at high values of delta T,
%% the results in this experiment are illustrated for two cases.
%% - One involving all 5000 data points.
%% - The second involving only the first 600 walks.
%% The experimental and theoretical values of diffusion coefficient are closer for the latter case.

N = 600;
MSD2 = zeros(N-1, 1);
deltaT2 = t(2:N) - t(1)*ones(N-1, 1);

% loop over the different possible time windows, calculate MSD for each
% time window size
for ii = 1:(N-1)
	
	%% -- Each loop deals with delT = 1, 2, 3, ... and computes one row of the MSD at the end.
	
	% calculate x displacements, separated by ii indices
	dx = x(1+ii:end) - x(1:end-ii);
	
	% calculate y displacements similarly
	dy = y(1+ii:end) - y(1:end-ii);
	
	% take mean over all displacements
	dispMean = mean( dx.^2 + dy.^2 );
	
	% store in MSD array
	MSD2(ii) = dispMean;
end
plot(deltaT2, MSD2, 'Color', 'green', 'LineWidth', 1.5, 'LineStyle', '--');
hold on;

%% (ADDITIONAL ends) ------------------------------------------------------------------------------


% get coefficients for line of first few MSD points using polyfit
linearCoeffs = polyfit(log10(deltaT), log10(MSD), 1);
% for first 600 points
linearCoeffs2 = polyfit(log10(deltaT2), log10(MSD2), 1);

% store slope and diffusion coefficient.
alpha = linearCoeffs(1);
diffusionCoefficient = 10^linearCoeffs(2);
% for first 600 points
alpha2 = linearCoeffs2(1);
diffusionCoefficient2 = 10^linearCoeffs2(2);

plot(deltaT, diffusionCoefficient*(deltaT.^alpha), 'Color', 'blue', 'LineWidth', 1.5, 'LineStyle', ':');
hold on;

plot(deltaT2, diffusionCoefficient2*(deltaT2.^alpha2), 'Color', 'black', 'LineWidth', 1.5, 'LineStyle', ':');
hold off;

legend('Experimental (all walks)', 'Experimental (first 600 walks)', 'Theoretical (all walks)', 'Theoretical (first 600 walks)');
fontsize(22, 'points');
xlabel('log(delT) - unitless');
ylabel('log(MSD) - unitless');

% print results to console
disp(['alpha is = ' num2str(alpha)]);
disp(['alpha is (based on first 600 walks) = ' num2str(alpha2)]);
disp(['measured diffusion coefficient is = ' num2str(diffusionCoefficient)]);
disp(['measured diffusion coefficient is (based on first 600 walks) = ' num2str(diffusionCoefficient2)]);


%% Part 5. Compare expected and measured diffusion coefficient

% store step size and step duration
stepSize = MSD(1);
stepDuration = t(2)-t(1);

% calculate expected diffusion coefficient
expectedD = stepSize*stepSize/stepDuration;

% print result to console
disp(['expected diffusion coefficient is = ' num2str(expectedD)]);
disp(" ");
disp('Since the amount of statistical error increases at high values of delta T, the results in this experiment are illustrated for two cases. One involving all 5000 data points. The second involving only the first 600 walks. The experimental and theoretical values of diffusion coefficient are closer for the latter case.');

%% PLEASE READ ADDITIONAL SECTION %%


