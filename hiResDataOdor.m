function [pulseIdx, pulseVol]=hiResDataOdor(dataFolder)

%% Navigate to 3dbrain folder and set as path

cd /projects/LEIFER/communalCode/3dbrain/
cd ..
path(pathdef)
cd 3dbrain/

% UI to pick BrainScanner folder

% if nargin==0
%     mostRecent=getappdata(0,'mostRecent');
% dataFolder=uipickfiles('FilterSpec',mostRecent);
% dataFolder=dataFolder{1};
% end

dataFolder = '/projects/LEIFER/PanNeuronal/20220516/BrainScanner20220516_160428';

%% Load other-flowdata.txt

%Load the other-flowdata file




