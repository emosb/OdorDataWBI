function []=pulse_behavior_overlay(dataFolder)

%% Navigate to 3dbrain folder and set as path

cd /projects/LEIFER/communalCode/3dbrain/
cd ..
path(pathdef)
cd 3dbrain/

% UI to pick BrainScanner folder

if nargin==0
    mostRecent=getappdata(0,'mostRecent');
dataFolder=uipickfiles('FilterSpec',mostRecent);
dataFolder=dataFolder{1};
setappdata(0,'mostRecent',fileparts(dataFolder));
end

% dataFolder = '/projects/LEIFER/PanNeuronal/20220527/BrainScanner20220527_113604';

%Load heatData

if exist([dataFolder filesep 'heatData.mat'],'file')
    heatData=load([dataFolder filesep 'heatData']);
    
end


%Store velocity
velocity = heatData.behavior.v;


%Load odorData

if exist([dataFolder filesep 'other-flowdata.txt'],'file')
    OdorData = readmatrix([dataFolder filesep 'other-flowdata.txt']);
    OdorDiff = diff(OdorData(:,3));
    pulseOn = find(OdorDiff>3);
    pulseOff = find(OdorDiff<-3);
    pulseVolOn = OdorData(pulseOn, 2);
    pulseVolOff = [1 OdorData(pulseOff,2)']';
    if pulseVolOn(end)>pulseVolOff(end)
        pulseVolOff = [pulseVolOff; length(velocity)];
    end
    
end

%Get pulse data sorted to MEK and M9
MEKpulsenum=length(pulseVolOn);
M9pulsenum=length(pulseVolOff)-1;
OdorAll = [];
for i=1:(max(MEKpulsenum,M9pulsenum))
    OdorAll = [OdorAll; pulseVolOff(i) pulseVolOn(i); pulseVolOn(i) pulseVolOff(i+1)];

end
M9pulseIdx = [];
MEKpulseIdx = [];
for i=1:(MEKpulsenum+M9pulsenum)
    if mod(i,2)==1
        M9pulseIdx = [M9pulseIdx; OdorAll(i,:)];
    else
        MEKpulseIdx = [MEKpulseIdx; OdorAll(i,:)];
    end
    
end

%Get rid of any indices beyond kept values in heatData
for i=1:MEKpulsenum
    if MEKpulseIdx(i,2)>length(velocity)
        MEKpulseIdx = MEKpulseIdx(1:i,:);
        MEKpulseIdx(i,2) = length(velocity);
        MEKpulsenum=length(MEKpulseIdx);
    end
end

for i=1:M9pulsenum
    if M9pulseIdx(i,2)>length(velocity)
        M9pulseIdx = M9pulseIdx(1:i,:);
        M9pulseIdx(i,2) = length(velocity);
        M9pulsenum=length(M9pulseIdx);
    end
end

%Fill in velocity values for MEK pulses
MEKPulseLength=MEKpulseIdx(:,2)-MEKpulseIdx(:,1);
MaxMEKPulseLength=max(MEKPulseLength);
vMEK=NaN(MaxMEKPulseLength+1, MEKpulsenum);

for i=1:MEKpulsenum
    vMEK(1:MEKPulseLength(i)+1,i)=velocity(MEKpulseIdx(i,1):MEKpulseIdx(i,2));
end

%Do the same for M9 pulses
M9PulseLength=M9pulseIdx(:,2)-M9pulseIdx(:,1);
MaxM9PulseLength=max(M9PulseLength);
vM9=NaN(MaxM9PulseLength+1, M9pulsenum);

for i=1:M9pulsenum
    vM9(1:M9PulseLength(i)+1,i)=velocity(M9pulseIdx(i,1):M9pulseIdx(i,2));
end

%Exclude first pulse - can exclude more or none by altering startPulse
startPulse=2;

%Find mean values for later plotting of mean velocity
vMeanM9=mean(vM9(:,startPulse:end),2,'omitnan');
vMeanMEK=mean(vMEK(:,startPulse:end),2,'omitnan');

%% Plot neural data and velocity with pulse location as vertical bars

%Limits for Plotting
pulseVol = [pulseVolOn; pulseVolOff];
y1=ones(length(pulseVol),1)*min(velocity);
y2=ones(length(pulseVol),1)*max(velocity);
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
imagesc(heatData.Ratio2)
%Plot vertical lines where flashes occur
xline(pulseVol, 'r')
xlim([pulseVol(1) length(heatData.hasPointsTime)])
title('Neural Data')

%Plotting velocity
nexttile(t)
plot(velocity)
hold on
%Plot vertical lines where flashes occur
plot([pulseVol, pulseVol]',[y1,y2]','r--')
%Match x-axis with neural data and only plot to max and min velocity to
%match vertical lines
xlim([pulseVol(1) length(heatData.hasPointsTime)])
ylim([min(velocity(pulseVol(1):end)) max(velocity(pulseVol(1):end))])
title('Velocity')

%Display dropped frames 
% disp(droppedFrames)

%% Plot the raw velocities

figure('Renderer', 'painters', 'Position', [1000 1000 1000 800])
t=tiledlayout(1,2,'TileSpacing','Compact');
title(t,'Odor Presence and Velocity')
xlabel(t,'Volume')
ylabel(t,'Velocity')

% M9 Pulse
nexttile(t)
hold on
for i=startPulse:M9pulsenum
    plot(vM9(:,i));
end
hold off
title('M9')

% MEK Pulse
nexttile(t)
hold on
for i=startPulse:MEKpulsenum
    plot(vMEK(:,i));
end
hold off
title('MEK')

%% Plotting average velocity over multiple pulses

figure('Renderer', 'painters', 'Position', [1000 1000 1000 800])
t2=tiledlayout(1,2,'TileSpacing','Compact');
title(t2,'Odor Presence and Mean Velocity')
xlabel(t2,'Volume')
ylabel(t2,'Mean Velocity Over First 45s of Cycle')
yMin=min(min(vMeanM9),min(vMeanMEK));
yMax=max(max(vMeanM9),max(vMeanMEK));

% M9 Pulse
nexttile(t2)
hold on
plot(vMeanM9);
yline(0, '--')
hold off
ylim([yMin yMax])
%Adjust number of volumes plotted (here 45 seconds roughly 270 volumes)
xlim([0 270])

title('M9')

% MEK Pulse
nexttile(t2)
hold on
plot(vMeanMEK);
yline(0, '--')
hold off
ylim([yMin yMax])
%Adjust number of volumes plotted
xlim([0 270])

title('MEK')

% %% Plotting only the first 30 seconds = roughly 180 volumes
% 
% figure('Renderer', 'painters', 'Position', [1000 1000 1000 800])
% t30=tiledlayout(2,2,'TileSpacing','Compact');
% title(t30,'Odor Presence and Velocity')
% xlabel(t30,'Volume')
% ylabel(t30,'Velocity')
% yMin=min(min(vMeanM9),min(vMeanMEK));
% yMax=max(max(vMeanM9),max(vMeanMEK));
% 
% % M9 Pulse Raw
% nexttile(t30)
% hold on
% for i=startPulse:M9pulsenum
%     plot(vM9(:,i));
% end
% yline(0,'--')
% hold off
% xlim([0, 180])
% title('M9 Raw')
% 
% % MEK Pulse Raw
% nexttile(t30)
% hold on
% for i=startPulse:MEKpulsenum
%     plot(vMEK(:,i));
% end
% yline(0,'--')
% hold off
% xlim([0, 180])
% title('MEK Raw')
% 
% % M9 Pulse Mean
% nexttile(t30)
% hold on
% plot(vMeanM9);
% yline(0, '--')
% hold off
% ylim([yMin yMax])
% xlim([0, 180])
% title('M9 Mean')
% 
% % MEK Pulse Mean
% nexttile(t30)
% hold on
% plot(vMeanMEK);
% yline(0, '--')
% hold off
% ylim([yMin yMax])
% xlim([0, 180])
% title('MEK Mean')

%% Calculate Neuronal Data of responsive neurons by pulse

sortedR2 = heatData.R2;
maxR2 = max(sortedR2, [], 2);
[~, index] = sort(maxR2, 'descend');
sortedR2 = sortedR2(index, :);

%Specify which neuron position in sorted list (e.g. 1 means most
%responsive, 2 means second most responsive, etc.)
neuron = 1;

neuronData=sortedR2(neuron,:)';

%Isolate data during M9 and MEK pulses
nM9=NaN(MaxM9PulseLength+1, M9pulsenum);

for i=1:M9pulsenum
    nM9(1:M9PulseLength(i)+1,i)=neuronData(M9pulseIdx(i,1):M9pulseIdx(i,2));
    
end

nMEK=NaN(MaxMEKPulseLength+1, MEKpulsenum);

for i=1:MEKpulsenum
    nMEK(1:MEKPulseLength(i)+1,i)=neuronData(MEKpulseIdx(i,1):MEKpulseIdx(i,2));
    
end

%Find mean values for later plotting of mean neuronal data
nMeanM9=mean(nM9(:,startPulse:end),2,'omitnan');
nMeanMEK=mean(nMEK(:,startPulse:end),2,'omitnan');

%% Plot Raw Data of specified neuron

figure('Renderer', 'painters', 'Position', [1000 1000 1000 800])
tnr=tiledlayout(1,2,'TileSpacing','Compact');
plotTitle=append("Odor Presence and Neuron ",num2str(index(neuron)),...
    " Response");
title(tnr,plotTitle)
xlabel(tnr,'Volume')
ylabel(tnr,'Neuronal Data')

startPulse=2;

% M9 Pulse
nexttile(tnr)
hold on
for i=startPulse:M9pulsenum
    plot(nM9(:,i));
end
hold off
title('M9')

% MEK Pulse
nexttile(tnr)
hold on
for i=startPulse:MEKpulsenum
    plot(nMEK(:,i));
end
hold off
title('MEK')

%% Mean Neuronal Data

figure('Renderer', 'painters', 'Position', [1000 1000 1000 800])
tnm=tiledlayout(1,2,'TileSpacing','Compact');
plotTitle=append("Odor Presence and Mean Neuron ",num2str(index(neuron)),...
    " Response Over First 45s of Cycle");
title(tnm,plotTitle)
xlabel(tnm,'Volume')
ylabel(tnm,'Mean Neuron Response')
yMin=min(min(nMeanM9),min(nMeanMEK));
yMax=max(max(nMeanM9),max(nMeanMEK));

% M9 Pulse
nexttile(tnm)
hold on
plot(nMeanM9);
yline(0, '--')
hold off
ylim([yMin yMax])
%Adjust number of volumes shown
xlim([0 270])

title('M9')

% MEK Pulse
nexttile(tnm)
hold on
plot(nMeanMEK);
yline(0, '--')
hold off
ylim([yMin yMax])
%Adjust number of volumes shown
xlim([0 270])

title('MEK')