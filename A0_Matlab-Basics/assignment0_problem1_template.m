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
nFlips = ;

% set number of trials
nTrials = ;

% use rand and round functions to generate each trial of coin flipping
% which will be represented by a nTrials x nFlips matrix of 1s and 0s


% calculate nHeads and nTails, the number of head and tail counts from each
% trial (hint: use sum function)
nHeads = ;
nTails = ;

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
figure(101),

hobj_heads = ;
hobj_tails = ;

% open figure window
figure(1), clf, hold on, box on;

% plot curves, label axes, add legend, etc.


%% Part 4. Plot multiple distributions as a function of nTrials

% repeat steps above to create three different distributions of head counts
% for different numbers of independent trials

% fix nFlips to be 100
nFlips = 100;



