function []=pulse_behavior_overlay_old(dataFolder)

%NOTE FOR NOW THIS SCRIPT ONLY WORKS FOR SPECIFIC DATA WHICH FOLLOWS
%SPECIFIC FLASHING PATTERNS - WILL ATTEMPT TO GENERALIZE LATER USING
%DISTANCE BETWEEN FLASH FRAMES > THRESHOLD

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

dataFolder = '/projects/LEIFER/PanNeuronal/20220516/BrainScanner20220516_151952';


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

%Store velocity
velocity = heatData.behavior.v;

%Split into ranges of butanone on and off
% OdorSwitch = volumeFlashLoc(3:16);
% OdorSwitch = [OdorSwitch(1:end-1)+1, OdorSwitch(2:end)-1];

duration = diff(volumeFlashLoc);
OdorPulse = find(duration>200);
OdorSwitch = [volumeFlashLoc(OdorPulse)+1, volumeFlashLoc(OdorPulse+1)-1];
OdorSwitch = OdorSwitch(find(OdorSwitch(:, 2)<length(velocity)),:);

MEKpulsenum=floor((length(OdorSwitch)-1)/2);
M9pulsenum=ceil(length(OdorSwitch)/2);
MEKpulseIdx = zeros(MEKpulsenum,2);
M9pulseIdx = zeros(M9pulsenum,2);
M9idx=1;
MEKidx=1;
for i=1:length(OdorSwitch)
    
    if mod(i,2)
        M9pulseIdx(M9idx,1) = OdorSwitch(i,1);
        M9pulseIdx(M9idx,2) = OdorSwitch(i,2);
        M9idx=M9idx+1;
    else
        MEKpulseIdx(MEKidx,1) = OdorSwitch(i,1);
        MEKpulseIdx(MEKidx,2) = OdorSwitch(i,2);
        MEKidx=MEKidx+1;
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
ylabel(t2,'Mean Velocity')
yMin=min(min(vMeanM9),min(vMeanMEK));
yMax=max(max(vMeanM9),max(vMeanMEK));

% M9 Pulse
nexttile(t2)
hold on
plot(vMeanM9);
yline(0, '--')
hold off
ylim([yMin yMax])
title('M9')

% MEK Pulse
nexttile(t2)
hold on
plot(vMeanMEK);
yline(0, '--')
hold off
ylim([yMin yMax])
title('MEK')

%% Plotting only the first 30 seconds = roughly 180 volumes

figure('Renderer', 'painters', 'Position', [1000 1000 1000 800])
t30=tiledlayout(2,2,'TileSpacing','Compact');
title(t30,'Odor Presence and Velocity')
xlabel(t30,'Volume')
ylabel(t30,'Velocity')
yMin=min(min(vMeanM9),min(vMeanMEK));
yMax=max(max(vMeanM9),max(vMeanMEK));

% M9 Pulse Raw
nexttile(t30)
hold on
for i=startPulse:M9pulsenum
    plot(vM9(:,i));
end
yline(0,'--')
hold off
xlim([0, 180])
title('M9 Raw')

% MEK Pulse Raw
nexttile(t30)
hold on
for i=startPulse:MEKpulsenum
    plot(vMEK(:,i));
end
yline(0,'--')
hold off
xlim([0, 180])
title('MEK Raw')

% M9 Pulse Mean
nexttile(t30)
hold on
plot(vMeanM9);
yline(0, '--')
hold off
ylim([yMin yMax])
xlim([0, 180])
title('M9 Mean')

% MEK Pulse Mean
nexttile(t30)
hold on
plot(vMeanMEK);
yline(0, '--')
hold off
ylim([yMin yMax])
xlim([0, 180])
title('MEK Mean')

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
    " Response");
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
title('M9')

% MEK Pulse
nexttile(tnm)
hold on
plot(nMeanMEK);
yline(0, '--')
hold off
ylim([yMin yMax])
title('MEK')



















