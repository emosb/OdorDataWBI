function []=pulse_alignment_old(dataFolder)

%Navigate to 3dbrain folder and set as path

cd /projects/LEIFER/communalCode/3dbrain/
cd ..
path(pathdef)
cd 3dbrain/

%UI to pick BrainScanner folder

if nargin==0
    mostRecent=getappdata(0,'mostRecent');
dataFolder=uipickfiles('FilterSpec',mostRecent);
dataFolder=dataFolder{1};
end

%Load hiResData

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
    
end

%Load heatData

if exist([dataFolder filesep 'heatData.mat'],'file')
    heatData=load([dataFolder filesep 'heatData']);

end

%Find volumes where high mag flash took place
volumeFlashLoc = hiResData.stackIdx(hiResData.flashLoc);

%Calculate number of dropped frames
droppedFrames = hiResData.stackIdx(end)-length(heatData.hasPointsTime);

%Store velocity
velocity = heatData.behavior.v;

%Limits for Plotting
y1=ones(length(volumeFlashLoc),1)*min(velocity);
y2=ones(length(volumeFlashLoc),1)*max(velocity);

%Alternative limits (to safely span vertical axis)
% y1=-100;
% y2=100;


%Plot the figure, neural data with flashes and velocity

figure('Renderer', 'painters', 'Position', [1000 1000 1600 1200])
t=tiledlayout(2,1,'TileSpacing','Compact');

%Grabbing plot title: date and file number

dataName=extractAfter(dataFolder,'BrainScanner');
dataString=convertCharsToStrings(dataName);
%Extract any strings of numbers, first is always date, second is file
%number
pat=digitsPattern;
dateAndFile=extract(dataString, pat);
plotTitle=append("Date:", dateAndFile(1), "  ","File:", dateAndFile(2));
title(t,plotTitle)
xlabel(t,'Volume Index')

% Plotting Neural Data
nexttile(t)
imagesc(heatData.R2)
%Plot vertical lines where flashes occur
xline(volumeFlashLoc, 'r')
xlim([volumeFlashLoc(1) length(heatData.hasPointsTime)])
title('Neural Data')

%Plotting velocity
nexttile(t)
plot(heatData.behavior.v)
hold on
%Plot vertical lines where flashes occur
plot([volumeFlashLoc, volumeFlashLoc]',[y1,y2]','r--')
%Match x-axis with neural data and only plot to max and min velocity to
%match vertical lines
xlim([volumeFlashLoc(1) length(heatData.hasPointsTime)])
ylim([min(velocity(volumeFlashLoc(1):end)) max(velocity(volumeFlashLoc(1):end))])
title('Velocity')

%Display dropped frames 
disp(droppedFrames)


