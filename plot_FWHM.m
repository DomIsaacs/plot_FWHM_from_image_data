clc;
close all;
clear all;

%% directory and analyis files

%file_n = matrix where rows are the pixel values along the line and columns
%represent individual frames. Typically these are stored as a CSV either
%from imageJ/FIJI directly or from cp into excel. 


%vessel_paths = ["C:/Users/......../animal/vessel/"];
%fileArray = ["file_1","file_2"....."file_n",];





%% imaging parameter
microns_per_pixel =(0.34); % microns/pixel
time_per_frame = 0.33; %s; 
StimFrameArray = [20]; % ---- --Where were the stimulations
number_of_baseline_frames = 19; % what is the "baseline"?



%% --------- diameter trace maker -----
for vessel_iter = 1:length(vessel_paths)
for file_iter = 1:length(fileArray)
filename = fileArray(file_iter);
path = vessel_paths(vessel_iter);
% ---------------------------
filepath = strcat(path,filename);
DATA = readtable(filepath);

DATA = table2array(DATA);

dataStart = 3; % column where data starts (its 3);
DATA = DATA(:,dataStart:end); %% scrap the first two columns no data when exporting the line profile from imageJ. 
Baseline_Diameter = [1:StimFrameArray(1)];
profile = DATA(:,[1 : end]);
End = size(profile,2);
detection_algorithm = 0;

%% using peak finder on background subtracted, time-grouped, edge-detected, smoothed images
% requires preprocessing (in ImageJ): 3x3 smoothing (Process->smooth),
% Rolling Ball Background Subtraction (~50pixel radius), Blurring Gaussian
% Kernel, then Edge Finder (alpha = 0.1-0.5) or just fwhm and extract line
% profile for every frame. Use the average intensity projection to pick the
% largest diameter 

if detection_algorithm == 1
ThisProfilePeakLocs = [];
junk = [];
x0 = profile;
for i = 1:size(x0,2)
    [ThisProfilePeakLocs, junk] = peakfinder(x0(:,i));
    if size(ThisProfilePeakLocs,1) == 2 % if it resolves two peaks indicating the Canny Edge Algorithm resolved a diameter 'edge'
         diameter(i) = microns_per_pixel*(ThisProfilePeakLocs(end,1) - ThisProfilePeakLocs(1,1));
    else
         diameter(i) = NaN;  % if it finds three peaks dont even try
    end
end
end


    
%% FWHM
if detection_algorithm == 0
    x = linspace(0,1*size(profile,1),size(profile,1));
    for i = 1:End
        xq = linspace(1,1*size(profile,1),1000);
        y = interp1(x',profile(:,i),xq','spline');
        FWHM(i) = fwhm(xq,y);   % <- MATLAB file exchange   
    end
    diameter = microns_per_pixel*FWHM;
end


%------------
%drop outlier frames (tends to be z-jitter or contrast dye uniformity disrupted)
std_FWHM_total = std(diameter);
TF = isoutlier(diameter,'mean'); % any data point higher than 3 std away from the mean of the whole time course. change to movemean
diameter(TF == 1) = NaN; % this method is poor becuase it will drop frames less at more variable (i.e. big responses that return to baseline) time series.
%------------------------


% -------------------------------------------------------------------------
% basal diameter (-x frames from 1st stim)
Do_span = StimFrameArray(1) - number_of_baseline_frames;
Do = mean(diameter(Do_span:StimFrameArray(1)),'omitnan');
sigma_o = std(diameter(Do_span:StimFrameArray(1)), 'omitnan'); % entire baseline variation.
% -------------------------------------------------------------------------
dataOut = diameter;

%diameter measurements in microns 
diameter_trace = diameter';
z_scored = (diameter - Do)/sigma_o;
percentChange = 100*(diameter - Do)/Do;
data.DiameterAcrossTime = diameter;
data.zscores = z_scored;
data.percentChange = percentChange;
data.Do = Do;
data.sigma_o = sigma_o;

data.baseline = z_scored(1:StimFrameArray(1))';

if length(StimFrameArray) == 1
data.S1 = z_scored(-number_of_baseline_frames + StimFrameArray(1)+1:end)';
end


save(strcat(path,strcat('OUTPUT_',filename,'.mat')),'data');
%% ------------------------FIRST SUBPLOT, FWHM, % change, z score to baseline mean and std

movmeanweight = 3;
figure;
subplot(3,1,1)
plot(data.DiameterAcrossTime,'k','linewidth',1);
M = movmean(data.DiameterAcrossTime,movmeanweight);
hold on
plot(M,'b','linewidth',2);
hold on
plot(Do.*ones(length(data.DiameterAcrossTime),1),'k--')
ylabel('microns');
xlim([-1, End*1.1]);
set(gca,'fontsize', 14)
title('profile FWHM');
xlabel('frame');
hold off


% 
subplot(3,1,2)
plot(percentChange,'k','linewidth',1);
M = movmean(percentChange,movmeanweight);
hold on
plot(M,'r','linewidth',2);
hold on
plot(zeros(length(data.DiameterAcrossTime),1),'k--')
ylabel('% change');
%ylim([-10 30]);
xlim([-1, End*1.1]);
set(gca,'fontsize', 14)
title('% change')
xlabel('frame');
hold off



subplot(3,1,3);
plot(z_scored,'k','linewidth',1);
M = movmean(z_scored,movmeanweight);
hold on
plot(M,'m','linewidth',2);
hold on
plot(zeros(length(data.DiameterAcrossTime),1),'k--');
ylabel('SD');
%ylim([-5 5]);
xlim([-1, End*1.1]);
set(gca,'fontsize', 14)
title('z-score') 
xlabel('frame');
hold off




%% ----------
% timeseries, calculations, etc.

movmeanweight_output = 1; %

allData.percentChangeTraces(:,file_iter) = movmean(percentChange',movmeanweight_output);
allData.DiameterTraces(:,file_iter) = movmean(diameter',movmeanweight_output);
allData.percentChangeRAW(:,file_iter) = percentChange';
allData.DiameterTracesRAW(:,file_iter) = diameter';
allData.DiameterDeltaRaw(:,file_iter) = diameter' - Do;
allData.FractionalChange(:,file_iter) = diameter/Do;
allData.zscores(:,file_iter) = movmean(z_scored',movmeanweight_output)';
allData.sd(:,file_iter) = sigma_o;
allData.Do(:,file_iter) = Do;



end

allData.percentChangeAveraged = mean(allData.percentChangeTraces,2);
allData.dilation = mean(allData.eachTrial_delta);
allData.baseline = mean(allData.Do);

end
