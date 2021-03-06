%% Description

% % Code for calculating the averaged metrics of MHWs by years for each location. Metrics includes frequency, 
%   total number of MHW days, mean duration, averaged mean intensity, mean maximum intensity,
%   mean cumulative intensity and mean intensity variability (standard deviation) 
%   of all MHWs in each year. In each location, the years with no MHWs have
%   the zero value, '0', for all columns


%% Input arguments:

%  start_yr - the first year of the study period [YYYY]
%  end_yr - the last year of the study period [YYYY]

%% Output arguments:

%  MHW_Annual_table - A table array containing the yearly averaged metrics
%  of MHWs where each row represents the annual averaged MHW metrics of one
%  location of a certain year. The averaged MHW metrics include:
%       - 'MHW_year' - one of the years during the study period [YYYY]
%       - 'Frequency' - total number of MHW events in one year [counts]
%       - 'Day' - total number of days with MHW in one year [days]
%       - 'MeanDuration' - mean duration of all MHWs in one year [days]
%       - 'MeanInt' - averaged mean intensity of all MHWs in one year [deg. C]
%       - 'MeanCumInt' - averaged cumulative intensity of all MHWs in one year [deg. C]
%       - 'MeanMaxInt' - mean maximum intensity of all MHWs in one year [deg. C]
%       - 'MeanVarInt' - averaged variance of intensity of all MHWs in one year [deg. C]
%       - 'Longitude' - longitude of the location
%       - 'Latitude' - latitude of the location

%  MHW_Annual - A double precision matrix array of MHW_Annual_table


%% Housekeeping
clc
clear all

%% Load array of all MHW events 
load MHW_All

%% Input needed
start_yr = 1983;
end_yr = 2012;

%% Find location of all MHWs
loc_full = unique(MHW_All(:,10:11),'rows');

%% Annual average of MHW metrics
MHW_Annual = [];

for n = 1:length(loc_full) % 0 to 90E
    [D1] = find(MHW_All(:,10) == loc_full(n,1) & MHW_All(:,11) == loc_full(n,2));
    for year = start_yr:end_yr
        [D2] = find(MHW_All(D1,8) == year);
        if ~isempty(D2)            
            Annual_Frequency = length(D2);
            Annual_Days = sum(MHW_All(D1(D2),3));
            Mean_Duration = mean(MHW_All(D1(D2),3),'all');
            Mean_MaxInt = mean(MHW_All(D1(D2),4),'all');
            Mean_MeanInt = mean(MHW_All(D1(D2),5),'all');
            Mean_VarInt = mean(MHW_All(D1(D2),6),'all');
            Mean_CumInt = mean(MHW_All(D1(D2),7),'all');                                                    
        else 
            Annual_Frequency = 0;
            Annual_Days = 0;
            Mean_Duration = 0;
            Mean_MeanInt = 0;
            Mean_CumInt = 0;
            Mean_MaxInt = 0;
            Mean_VarInt = 0;                                   
        end
        
        MHW_year = year;
        longitude = loc_full(n,1);
        latitude = loc_full(n,2);
            
        MHW_Annual = [MHW_Annual;[MHW_year Annual_Frequency Annual_Days Mean_Duration Mean_MeanInt Mean_CumInt Mean_MaxInt Mean_VarInt longitude latitude]];
    end    
end

%% Create table array from matrix
MHW_Annual_table = array2table(MHW_Annual, 'VariableNames', {'MHW_year', 'Frequency', 'Day', 'MeanDuration', 'MeanInt', 'MeanCumInt', 'MeanMaxInt', 'MeanVarInt', 'Longitude','Latitude'});

%% Save MHW_Annual to directory path for further analysis
save("MHW_Annual","MHW_Annual");

%% %% Clear un-needed variables
clearvars -except MHW_All MHW_Annual MHW_Annual_table

%% Concatenate MHW_Annual if needed
% load MHW_Annual1 MHW_Annual2 MHW_Annual3      %% Example filename
% MHW_All = [MHW_Annual1;MHW_Annual2;MHW_Annual3];
