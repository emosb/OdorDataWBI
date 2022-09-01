function []=reversal_time_by_odor(dataFolder)

%% Navigate to the BrainScanner folder

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

% Load heatData

if exist([dataFolder filesep 'heatData.mat'],'file')
    heatData=load([dataFolder filesep 'heatData']);
    
end


% Store velocity
velocity = heatData.behavior.v;


% Load odorData

if exist([dataFolder filesep 'other-flowdata.txt'],'file')
    OdorData = readmatrix([dataFolder filesep 'other-flowdata.txt']);
    pulseOn = OdorData(find(diff(OdorData(:,3))>3),2);
    pulseOff = OdorData(find(diff(OdorData(:,3))<-3),2);
    
end

% Sort pulse volumes into a full list with labels
OdorSwitch = [pulseOn zeros(size(pulseOn)); pulseOff ones(size(pulseOff))];
switchVol = zeros((length(pulseOn)+length(pulseOff)),3);
switchVol(1,1)=1;
for i=1:(length(pulseOn)+length(pulseOff))
    [pulseEnd, Idx] = min(OdorSwitch(:,1));
    if pulseEnd>=length(velocity)
        switchVol(end,2) = length(velocity);
        break;
    end
    switchVol(i,2:end) = OdorSwitch(Idx,:);
    switchVol(i+1,1) = pulseEnd+1; 
    OdorSwitch(Idx,1) = NaN;
    if i==(length(pulseOn)+length(pulseOff))
        switchVol(end,2) = length(velocity);
    end
end

%Exclude first pulse - can exclude more or none by altering startPulse
startPulse=2;

%Find mean and stdev of all velocities excluding omitted pulse(s)
startVol = switchVol((2*startPulse)-1,1);
endVol = length(velocity);
vMeanAll = mean(velocity(startVol:endVol),'omitnan');
vStdAll = std(velocity(startVol:endVol),'omitnan');

%We want a moving threshold for each full cycle - to account for changing
%stddev of velocity as recording progresses
if switchVol(((2*startPulse)-1),3) == switchVol(end,3)
    endIdx = length(switchVol(((2*startPulse)-1):end-1,1));
else
    endIdx = length(switchVol(((2*startPulse)-1):end),1);
end 

%for each individual cycle, calculate mean and std of velocity
vMeanCycle = zeros((floor(endIdx/2)),1);
vStdCycle = zeros((floor(endIdx/2)),1);
vThresCycle = zeros((floor(endIdx/2)),1);
RevPoints = zeros((floor(endIdx/2)),2);
j = (2*startPulse)-1;
for i=1:(floor(endIdx/2))
    vMeanCycle(i) = mean(velocity(switchVol(j,1):switchVol(j+1,2)), 'omitnan');
    vStdCycle(i) = std(velocity(switchVol(j,1):switchVol(j+1,2)), 'omitnan');
    vThresCycle(i) = vMeanCycle(i) - vStdCycle(i);
    for k=switchVol(j,1):switchVol(j,2)
       if velocity(k)<vThresCycle(i) 
           RevPoints(i,1) = RevPoints(i,1)+1;
       end
    end
    for k=switchVol(j+1,1):switchVol(j+1,2)
       if velocity(k)<vThresCycle(i) 
           RevPoints(i,2) = RevPoints(i,2)+1;
       end
    end
    j = j+2;
end

RevPointsMean = [mean(RevPoints(:,1), 'omitnan') mean(RevPoints(:,2), 'omitnan')];
RevPointsStd = [std(RevPoints(:,1), 'omitnan') std(RevPoints(:,2), 'omitnan')];

%Bar plot with error bars
x=1:2;
errhigh = RevPointsStd;
errlow  = RevPointsStd;

bar(x,RevPointsMean)   
set(gca,'xticklabel',{'M9'; 'MEK'})

hold on

er = errorbar(x,RevPointsMean,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

% %Create threshold of velocities below X standard deviations from the mean
% vThreshold = vMean - 1*vStd; %Possible we replace stddev with some other metric later
% RevPoint = zeros(length(velocity),2); %Keep track of reversal volume indices
% for i=startVol:endVol
%    if velocity(i)<vThreshold
%        RevPoint(i,1) = 1;
%    end
% end

%Plot the reversal velocities (just for visualization also)
% figure;
% plot(velocity.*RevPoint)

% %Bin reversal volumes into odor type (M9 or MEK)
% RevPoint=find(RevPoint(:,1)==1);
% revTimeM9 = zeros(length(switchVol),1);
% revTimeMEK = zeros(length(switchVol),1);
% for i=1:length(switchVol)
%     if i<=startPulse
%         revTimeM9(i) = NaN;
%         revTimeMEK(i) = NaN;
%     else
%         for j=1:length(RevPoint)
%             if (RevPoint(j)>switchVol(i,1) && RevPoint(j)<switchVol(i,2))
%                 if switchVol(i,3)==1
%                     revTimeMEK(i) = revTimeMEK(i) +1;
%                     revTimeM9(i) = NaN;
%                 else
%                     revTimeM9(i) = revTimeM9(i) +1;
%                     revTimeMEK(i) = NaN;
%                 end
%             end
%         end
%     end
% end
% 
% meanRevTimeMEK = mean(revTimeMEK, 'omitnan')/length(RevPoint);
% stdRevTimeMEK = std(revTimeMEK, 'omitnan')/length(RevPoint);
% meanRevTimeM9 = mean(revTimeM9, 'omitnan')/length(RevPoint);
% stdRevTimeM9 = std(revTimeM9, 'omitnan')/length(RevPoint);
% 
% % figure;
% % bar([revTimeM9/length(RevPoint) revTimeMEK/length(RevPoint)])
% % set(gca,'xticklabel',{'M9'; 'MEK'})
% 
% x = 1:2;
% data = [meanRevTimeM9; meanRevTimeMEK];
% errhigh = [stdRevTimeM9; stdRevTimeMEK];
% errlow  = [stdRevTimeM9; stdRevTimeMEK];
% 
% bar(x,data)   
% set(gca,'xticklabel',{'M9'; 'MEK'})
% 
% hold on
% 
% er = errorbar(x,data,errlow,errhigh);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% hold off

