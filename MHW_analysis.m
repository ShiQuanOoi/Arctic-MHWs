%% Description

% % Code for analysing all the detected MHWs. Analysis is done on MHW
%   frequency, total days, mean duration and mean intensity only. Different
%   types of analysis are included and are seperated into sections:
%       1. MHW averaged metrics over 30 years - Calculate averaged MHW frequency, 
%          total days, duration and mean intensity over 30 years

%       2. Linear trends of MHW metrics by year - Calculate linear trend of MHW frequency, 
%          total days, duration and mean intensity by year

%       3. Differences of MHW metrics between 1983-1997 and 1998-2012
%          periods - Calculate differences of MHW mean frequency, mean total days, mean
%          duration and mean intensity between 1983-1997 and 1998-2012 periods

%       4. Area-weighted averaged time series of annual MHW metrics - Plot
%          the time series of annual MHW mean frequency, mean total days, mean
%          duration and mean intensity in 3 latitude bands (60N-70N, 70N-80N and 80N-90N)

%       5. Histogram of frequency of each MHW duration between 1983-1997
%          and 1998-2012 periods - Compare the occurrence of each MHW
%          duration between 1983-1997 and 1998-2012 periods

%       6. Survivorship analysis on the frequency of each MHW duration
%          between 1983-1997 and 1998-2012 periods - Compare the
%          termination rate of the MHW events with different durations 
%          between 1983-1997 and 1998-2012 periods

%       7. Frequency of long MHWs (duration>14days) - Number of MHW events
%          with duration more than 14 days. This analysis consists of total
%          number of long MHWs over 30 years and number of long MHWs before
%          and after 1998.
  

% % Output arguments:

%   MHW_30yearsMean_table - Averaged MHW frequency, total days, duration
%   and mean intensity over 30 years

%   MHW_Trend_table - Linear trend of MHW frequency, total days, duration and mean
%   intensity by year

%   MHW_MeanDiff_table - Differences of MHW mean frequency, mean total days, mean
%   duration and mean intensity between 1983-1997 and 1998-2012 periods

%   MHW_14dFreq_table - Total number of long MHWs over 30 years

%   MHW_14dPre98Freq_table - Total number of long MHWs during 1983-1997

%   MHW_14dPost98Freq_table - Total number of long MHWs during 1998-2012



% % An additional function, 'trend', written by Chad Greene, is needed for calculation and can be
%   downloaded from https://uk.mathworks.com/matlabcentral/fileexchange/46363-trend


%% Housekeeping
clc
clear all
close all

%% Load array of all MHW events 
load MHW_All.mat
load MHW_Annual.mat

%% Find location of all MHWs
loc_full = unique(MHW_All(:,10:11),'rows');

%% Selects 4 metrics for analysis (frequency, total days, duration and mean intensity)
MHW_metrics = MHW_Annual(:,2:5);





%% 1. MHW averaged metrics over 30 years
MHW_30yearsMean = [];

for col_no = 1:size(MHW_metrics,2)
    for n = 0:length(loc_full)-1
        D1 = 1+n*30;
        D2 = 30+n*30;
        MHW_30yearsMean(n+1,col_no) = mean(MHW_metrics(D1:D2,col_no));
    end    
end

% Concatenate with location
MHW_30yearsMean = [MHW_30yearsMean loc_full];   

% Create table array from matrix
MHW_30yearsMean_table = array2table(MHW_30yearsMean,'VariableNames', {'Average_Frequency',...
    'Average_Day','Average_MeanDuration','Average_MeanInt','Longitude','Latitude'});

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table loc_full MHW_metrics






%% 2. Linear trends of MHW metrics by year
MHW_Trend = [];

for col_no = 1:size(MHW_metrics,2)
    for n = 0:length(loc_full)-1
        D1 = 1+n*30;
        D2 = 30+n*30;
        MHW_Trend(n+1,col_no) = trend(MHW_metrics(D1:D2,col_no));
    end    
end

% Concatenate with location
MHW_Trend = [MHW_Trend loc_full];

% Create table array from matrix
MHW_Trend_table = array2table(MHW_Trend,'VariableNames', {'T_Frequency','T_Day',...
    'T_MeanDuration','T_MeanInt','Longitude','Latitude'});

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table MHW_Trend_table loc_full 






%% 3. Differences of MHW metrics between 1983-1997 and 1998-2012 periods

% MHW averaged metrics in 1983-1997
MHW_AnnualPre98 = MHW_Annual(MHW_Annual(:,1)<1998,:);
MHW_metricsPre98 = MHW_AnnualPre98(:,2:5);

MHW_MeanPre98 = [];

for col_no = 1:size(MHW_metricsPre98,2)
    for n = 0:length(loc_full)-1
        D1 = 1+n*15;
        D2 = 15+n*15;
        MHW_MeanPre98(n+1,col_no) = mean(MHW_metricsPre98(D1:D2,col_no));
    end    
end


% MHW averaged metrics in 1998-2012
MHW_AnnualPost98 = MHW_Annual(MHW_Annual(:,1)>1997,:);
MHW_metricsPost98 = MHW_AnnualPost98(:,2:5);

MHW_MeanPost98 = [];

for col_no = 1:size(MHW_metricsPost98,2)
    for n = 0:length(loc_full)-1
        D1 = 1+n*15;
        D2 = 15+n*15;
        MHW_MeanPost98(n+1,col_no) = mean(MHW_metricsPost98(D1:D2,col_no));
    end    
end


% Calculate differences
MHW_MeanDiff = MHW_MeanPost98 - MHW_MeanPre98;
MHW_MeanDiff = [MHW_MeanDiff loc_full];
MHW_MeanDiff_table = array2table(MHW_MeanDiff,'VariableNames', {'Diff_Frequency','Diff_Day',...
    'Diff_MeanDuration','Diff_MeanInt','Longitude','Latitude'});

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table MHW_Trend_table MHW_MeanDiff_table 






%% 4. Area-weighted averaged time series of annual MHW metrics

% Seperate MHW_Annual array into 3 latitude bands (low: 60°N-70°N, middle: 70°N-80°N and high: 80°N-90°N)
MHW_annual_60N = MHW_Annual(MHW_Annual(:,10)<70,:);
MHW_annual_70N = MHW_Annual(MHW_Annual(:,10)<80 & MHW_Annual(:,10)>70,:);
MHW_annual_80N = MHW_Annual(MHW_Annual(:,10)>=80,:);

% Location array of 3 latitude bands
loc_60N = unique(MHW_annual_60N(:,9:10),'rows');
loc_70N = unique(MHW_annual_70N(:,9:10),'rows');
loc_80N = unique(MHW_annual_80N(:,9:10),'rows');

% Array with all years in the study period
year_full = MHW_Annual(1:30,1);

% Time series arrays
MHW_TimeSeries_60N = [];
MHW_TimeSeries_70N = [];
MHW_TimeSeries_80N = [];

MHW_TimeSeries_60N(:,1) = year_full;
MHW_TimeSeries_70N(:,1) = year_full;
MHW_TimeSeries_80N(:,1) = year_full;

for n = year_full(1):year_full(end)
    
    % averaged frequency
    MHW_TimeSeries_60N(n-1982,2) = mean(MHW_annual_60N(MHW_annual_60N(:,1)==n,2));
    MHW_TimeSeries_70N(n-1982,2) = mean(MHW_annual_70N(MHW_annual_70N(:,1)==n,2));
    MHW_TimeSeries_80N(n-1982,2) = mean(MHW_annual_80N(MHW_annual_80N(:,1)==n,2));
    
    % averaged total days
    MHW_TimeSeries_60N(n-1982,3) = mean(MHW_annual_60N(MHW_annual_60N(:,1)==n,3));
    MHW_TimeSeries_70N(n-1982,3) = mean(MHW_annual_70N(MHW_annual_70N(:,1)==n,3));
    MHW_TimeSeries_80N(n-1982,3) = mean(MHW_annual_80N(MHW_annual_80N(:,1)==n,3));
    
    % averaged duration
    MHW_TimeSeries_60N(n-1982,4) = mean(MHW_annual_60N(MHW_annual_60N(:,1)==n,4));
    MHW_TimeSeries_70N(n-1982,4) = mean(MHW_annual_70N(MHW_annual_70N(:,1)==n,4));
    MHW_TimeSeries_80N(n-1982,4) = mean(MHW_annual_80N(MHW_annual_80N(:,1)==n,4));
    
    % averaged mean intensity
    MHW_TimeSeries_60N(n-1982,5) = mean(MHW_annual_60N(MHW_annual_60N(:,1)==n,5));
    MHW_TimeSeries_70N(n-1982,5) = mean(MHW_annual_70N(MHW_annual_70N(:,1)==n,5));
    MHW_TimeSeries_80N(n-1982,5) = mean(MHW_annual_80N(MHW_annual_80N(:,1)==n,5));
    
end

% Create table arrays from matrices
MHW_TimeSeries_60N_T = array2table(MHW_TimeSeries_60N,'VariableNames', {'Year','MeanFrequency',...
    'MeanDay','MeanDuration','MeanInt'});
MHW_TimeSeries_70N_T = array2table(MHW_TimeSeries_70N,'VariableNames', {'Year','MeanFrequency',...
    'MeanDay','MeanDuration','MeanInt'});
MHW_TimeSeries_80N_T = array2table(MHW_TimeSeries_80N,'VariableNames', {'Year','MeanFrequency',...
    'MeanDay','MeanDuration','MeanInt'});


% Plotting time series
% Averaged frequency
plot(MHW_TimeSeries_60N(:,1),MHW_TimeSeries_60N(:,2),'k','LineWidth',1.1)
xlim([1983 2012])
xlabel('Year')
ylabel('Frequency (counts)')
hold on
plot(MHW_TimeSeries_70N(:,1),MHW_TimeSeries_70N(:,2),'b','LineWidth',1.1)
plot(MHW_TimeSeries_80N(:,1),MHW_TimeSeries_80N(:,2),'r','LineWidth',1.1)
hold off
legend('low','middle','high','FontSize',9)
figure

% Averaged total days
plot(MHW_TimeSeries_60N(:,1),MHW_TimeSeries_60N(:,3),'k','LineWidth',1.1)
xlim([1983 2012])
xlabel('Year')
ylabel('Total days (days)')
hold on
plot(MHW_TimeSeries_70N(:,1),MHW_TimeSeries_70N(:,3),'b','LineWidth',1.1)
plot(MHW_TimeSeries_80N(:,1),MHW_TimeSeries_80N(:,3),'r','LineWidth',1.1)
hold off
legend('low','middle','high','FontSize',9)
figure

% Averaged duration
plot(MHW_TimeSeries_60N(:,1),MHW_TimeSeries_60N(:,4),'k','LineWidth',1.1)
xlim([1983 2012])
xlabel('Year')
ylabel('Duration (days)')
hold on
plot(MHW_TimeSeries_70N(:,1),MHW_TimeSeries_70N(:,4),'b','LineWidth',1.1)
plot(MHW_TimeSeries_80N(:,1),MHW_TimeSeries_80N(:,4),'r','LineWidth',1.1)
hold off
legend('low','middle','high','FontSize',9)
figure

% Averaged mean intensity
plot(MHW_TimeSeries_60N(:,1),MHW_TimeSeries_60N(:,5),'k','LineWidth',1.1)
xlim([1983 2012])
xlabel('Year')
ylabel('Mean intensity ({\circ}C)')
hold on
plot(MHW_TimeSeries_70N(:,1),MHW_TimeSeries_70N(:,5),'b','LineWidth',1.1)
plot(MHW_TimeSeries_80N(:,1),MHW_TimeSeries_80N(:,5),'r','LineWidth',1.1)
hold off
legend('low','middle','high','FontSize',9)
figure

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table MHW_Trend_table MHW_MeanDiff_table
    






%% 5. Histogram of frequency of each MHW duration between 1983-1997 and 1998-2012 periods

% Find durations of all MHWs
MHW_duration = MHW_All(:,3);
duration_full = unique(MHW_duration);

% Seperate all MHW events into two periods, 1983-1997 and 1998-2012
MHW_AllPre98 = MHW_All(MHW_All(:,8)<1998,:);
MHW_AllPost98 = MHW_All(MHW_All(:,8)>=1998,:);

% Calculate the occurrence of each duration during pre and post 1998
FrequencyOfDuration_pre98 = [duration_full,histc(MHW_AllPre98(:,3),duration_full)];
FrequencyOfDuration_post98 = [duration_full,histc(MHW_AllPost98(:,3),duration_full)];


% Plot bar chart
bar(duration_full,FrequencyOfDuration_post98(:,2))
set(gca,'YScale','log')
hold on 
bar(duration_full,FrequencyOfDuration_pre98(:,2),'FaceAlpha',0.7)
xlabel('Duration (days)')
ylabel('Frequency (log10)')
legend('Post-1998','Pre-1998')
hold off
figure

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table MHW_Trend_table MHW_MeanDiff_table ...
    FrequencyOfDuration_pre98 FrequencyOfDuration_post98 duration_full






%% 6. Survivorship analysis on the frequency of each MHW duration between 1983-1997 and 1998-2012 periods

% Calculate cumulative frequency of the durations for two periods, 1983-1997 and 1998-2012
CumFreq_pre98 = cumsum(FrequencyOfDuration_pre98(:,2))/sum(FrequencyOfDuration_pre98(:,2));
CumFreq_post98 = cumsum(FrequencyOfDuration_post98(:,2))/sum(FrequencyOfDuration_post98(:,2));

% Plot
yy1 = 1-CumFreq_pre98;
yy2 = 1-CumFreq_post98;
semilogy(duration_full,yy1);
hold on
semilogy(duration_full,yy2);
hold off
xlabel('Duration (days)')
ylabel('log(1-P)')
legend('Pre-1998','Post-1998')

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table MHW_Trend_table MHW_MeanDiff_table 







%% 7. Frequency of long MHWs (duration>14days)

% Select long MHWs from all MHW events
MHW_14d = MHW_All(MHW_All(:,3)>14,:);
MHW_14d = MHW_14d(:,[3 8 10 11]);
loc_full = unique(MHW_14d(:,3:4),'rows');

% Calculate the total frequency of long MHWs over 30 years
MHW_14dFreq = [];

for n = 1:length(loc_full)
    [D1] = find(MHW_14d(:,3) == loc_full(n,1) & MHW_14d(:,4) == loc_full(n,2));
    MHW_14dFreq(n,:) = [length(D1) loc_full(n,1) loc_full(n,2)];
end
    
% Create table array from matrix
MHW_14dFreq_table = array2table(MHW_14dFreq,'VariableNames', {'Frequency','Longitude','Latitude'});


% Calculate the total frequency of long MHWs for two periods, 1983-1997 and 1998-2012
MHW_14dPre98 = MHW_14d(MHW_14d(:,2)<1998,:);
MHW_14dPost98 = MHW_14d(MHW_14d(:,2)>1997,:);

MHW_14dPre98Freq = [];
MHW_14dPost98Freq = [];

for n = 1:length(loc_full)
    [D1] = find(MHW_14dPre98(:,3) == loc_full(n,1) & MHW_14dPre98(:,4) == loc_full(n,2));
    if ~isempty(D1)
        MHW_14dPre98Freq(n,:) = [length(D1) loc_full(n,1) loc_full(n,2)];
    else 
        MHW_14dPre98Freq(n,:) = [0 loc_full(n,1) loc_full(n,2)];
    end
end

for n = 1:length(loc_full)
    [D1] = find(MHW_14dPost98(:,3) == loc_full(n,1) & MHW_14dPost98(:,4) == loc_full(n,2));
    if ~isempty(D1)
        MHW_14dPost98Freq(n,:) = [length(D1) loc_full(n,1) loc_full(n,2)];
    else 
        MHW_14dPost98Freq(n,:) = [0 loc_full(n,1) loc_full(n,2)];
    end
end

% Create table array from matrix
MHW_14dPre98Freq_table = array2table(MHW_14dPre98Freq,'VariableNames', {'Frequency','Longitude','Latitude'});
MHW_14dPost98Freq_table = array2table(MHW_14dPost98Freq,'VariableNames', {'Frequency','Longitude','Latitude'});

% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_30yearsMean_table MHW_Trend_table MHW_MeanDiff_table ...
    MHW_14dFreq_table MHW_14dPre98Freq_table MHW_14dPost98Freq_table

