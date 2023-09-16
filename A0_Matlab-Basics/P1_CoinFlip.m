fprintf("Script begins...\n");

%% Assignment 0, Problem 1 template

%{
    Problem 1. Coin Flipping Experiment

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

%% Part 1,2. Obtain number of heads and tails found after multiple trials of flipping multiple coins

% set number of coin flips
nFlips = 1000;

% set number of trials
nTrials = 100;

% use rand and round functions to generate each trial of coin flipping
% which will be represented by a nTrials x nFlips matrix of 1s and 0s
coinFlips = round(rand(nTrials, nFlips));

% calculate nHeads and nTails, the number of head and tail counts from each
% trial (hint: use sum function)

% assume: 1 is heads, 0 is tails.
nHeads = sum(coinFlips, 2);
nTails = nFlips*ones(nTrials, 1) - nHeads;


%% Part 3. Plot distributions of heads and tail counts

% =================================
% create distributions of the
% number of heads and tails
% found during trials
% =================================

% number of bins for histogram (can change as necessary)
nbins = 15;

% create histogram object using histogram function for head and tail flips, obtain distributions
% NOTE: histogram function will call a window, so run figure(101) before
% calling function to put output onto dummy window.
figure(101);

hobj_heads = histogram(nHeads,'Normalization','pdf');
heads_xcenters = hobj_heads.BinEdges(1:(end-1)) + hobj_heads.BinWidth/2;
heads_ys = hobj_heads.Values;

hobj_tails = histogram(nTails,'Normalization','pdf');
tails_xcenters = hobj_tails.BinEdges(1:(end-1)) + hobj_tails.BinWidth/2;
tails_ys = hobj_tails.Values;


figure(1);
hold on;
plot(heads_xcenters, heads_ys, "Color", 'blue', "LineWidth", 1.5);
hold on;
plot(tails_xcenters, tails_ys, "Color", "green", "LineWidth", 1.5);
hold off;
legend('heads', 'tails');
fontsize(22, 'points');
xlabel('N_{heads}');
ylabel('P(N_{heads})');

% open figure window
% figure(1), clf, hold on, box on;

% plot curves, label axes, add legend, etc.


%% Part 4. Plot multiple distributions as a function of nTrials

% repeat steps above to create three different distributions of head counts
% for different numbers of independent trials

% fix nFlips to be 100
nFlips = 100;
nTrialsVec = [100; 1000; 10000];
nColorsVec = ["blue"; "green"; "yellow"];

figure(2);
hold on;
for nTrialIdx = 1:size(nTrialsVec)
    
    nTrials = nTrialsVec(nTrialIdx);
    
    coinFlips = round(rand(nTrials, nFlips));
    nHeads = sum(coinFlips, 2);
    
    figure(101);
    hobj_heads = histogram(nHeads, 'Normalization', 'pdf');
    heads_xcenters = hobj_heads.BinEdges(1:(end-1)) + hobj_heads.BinWidth/2;
    heads_ys = hobj_heads.Values;
    
    figure(2);
    hold on;
    plot(heads_xcenters, heads_ys, "Color", nColorsVec(nTrialIdx), "LineWidth", 1.5);
    
end

figure(2);
hold on;
idealBinom_xs = 1:nFlips;
idealBinom_ys = binopdf(idealBinom_xs, nFlips, 0.5);
plot(idealBinom_xs, idealBinom_ys, "Color", "red", "LineWidth", 1.5, "LineStyle","--");
hold off;

plot()