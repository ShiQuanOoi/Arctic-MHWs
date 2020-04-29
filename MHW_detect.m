%% Description

% % The code for detecting marine heatwaves (MHWs), including loading SST data, calculate
%   climatological mean and threshold, recording every MHW in each
%   location and calculate the duration, mean intensity, maximum intensity,
%   cumulative intensity and intensity variability (standard deviation)
%   of each MHW event. 

% % The code is based on the MHW definition proposed by Hobday et al. (2016) 'A hierarchical 
    approach to defining marine heatwaves', Progress in Oceanography, 141, p227-238.

% % An additional function, 'findND', written by Rik, is needed for calculation and can be
%   downloaded from https://uk.mathworks.com/matlabcentral/fileexchange/64383-findnd


%% Input arguments:

%  YYYYmmdd - Start and end dates of the period (YYYY,MM,DD) 

%  YYYY - Start and end years of the period (YYYY,MM,DD)

%  lon_start - Starting point of the area in longitude (0.125 to 359.875 degree)

%  lat_start - Starting point of the area in latitude (-89.875 to 89.875 degree)

%  nLon - Define the size of the study area, number of grids needed in longitude

%  nLat - Define the size of the study area, number of grids needed in latitude


%% Output arguments:

%  MHW_table - A table containing all detected MHWs where each row
%  represents a single event and each column corresponds to a
%  metric. MHW metrics are: 
%       - 'MHW_Start' - start date of each event [YYYYMMDD]
%       - 'MHW_End' - end date of each event [YYYYMMDD]
%       - 'Duration' - duration of each event [days]
%       - 'MaxInt' - maximum intensity of each event [deg. C]
%       - 'MeanInt' - mean intensity of each event [deg. C]
%       - 'VarInt' - variance of intensity during each event [deg. C]
%       - 'CumInt' - cumulative intensity across each event [deg. C]
%       - 'MHW_StartYr' - starting year of each event [YYYY]
%       - 'MHW_EndYr' - ending year of each event [YYYY]
%       - 'Longitude' - longitude of each event 
%       - 'Latitude' - latitude of each event   

%  MHW_All - A double precision matrix array of MHW_table

%  SST_All = All SST data from the satellite record for the whole area and
%             whole study period in 4D matrix (lon by lat by 366 by year)

%  SST_Mean = Climatological mean of SST_All in 3D matrix (lon by lat by 366)

%  SST_90 = Climatological 90th percentile threshold of SST_All in 3D matrix (lon by lat by 366)


% % Un-needed variables are deleted regularly to save memory in the
%   workspace. Variables can be kept simply by disable the clearvars
%   function



%% Housekeeping
clc
clear all

%% Import parameters from OISST data, eg. location and time
% Read info from one of the dataset
f = ncinfo('sst.day.mean.1983.nc');
nvars = length(f.Variables);
for k = 1:nvars
    varname = f.Variables(k).Name;
    disp(['Reading:  ' varname]);
    eval([varname ' = ncread(''' 'sst.day.mean.1983.nc' ''',''' varname ''');']);
end

% Land value/ missing value
LandValue = f.Variables.FillValue;
LandValue = -LandValue;

% Global longitude and latitude arrays
latitudeAll = double(ncread('sst.day.mean.1983.nc','lat'));
longitudeAll = double(ncread('sst.day.mean.1983.nc','lon'));

%% Input needed
% Time array, put in start and end dates
YYYYmmdd = str2double(string(datestr(datenum(1983,1,1):datenum(2012,12,31),'YYYYmmdd')));   % format (yyyy,mm,dd)
YYYY = str2double(string(datestr(datenum(1983,1,1):datenum(2012,12,31),10)));

% Specify study area
% Find starting point from the longitude and latitude arrays
[lon_start,~] = find(longitudeAll==0.125);
[lat_start,~] = find(latitudeAll==60.125); 

% Define area by putting in number of longitude and latitude cells needed
nLon = 5;     
nLat = 5;    

%% Build seperate arrays represent all the longitude and latitude in the study area
lon_used = longitudeAll(lon_start:lon_start+nLon-1);
lat_used = latitudeAll(lat_start:lat_start+nLat-1);

%% Load SST data
% Define starting point and number of SST data needed for each year
start = [lon_start lat_start 1];    % 1 = first day of the year
count_nm = [nLon nLat 365];         % nm = normal year
count_lp = [nLon nLat 366];         % lp = leap year

% Set up a NaN row to insert in between 28 Feb and 1 Mar of normal year
leap_row = NaN(nLon,nLat,1);

% Set up array to store imported data
sst_extended = [];

% Import and concatenate
for i = 1982:2013     % Normal 365 days
    sst_nm = squeeze(ncread(sprintf('sst.day.mean.%d.nc',i),'sst',start,count_nm));
    sst_nm = cat(3,sst_nm(:,:,1:59),leap_row,sst_nm(:,:,60:end));   % Insert the leap row
    sst_extended(:,:,:,i-1981) = sst_nm; 
end

for k = 1984:4:2012   % Leap 366 days
    sst_extended(:,:,:,k-1981) = squeeze(ncread(sprintf('sst.day.mean.%d.nc',k),'sst',start,count_lp));
end

%% Replace land value/ missing value with NaN
sst_extended(sst_extended == LandValue) = NaN;

%% Clear un-needed variables
clearvars -except YYYY YYYYmmdd nLon nLat lat_used lon_used sst_extended

%% Make 11-day window average
% Insert first 5 days from the year later to the end of each year
for i = 3:32
    sst_long1(:,:,:,i-2) = cat(3,sst_extended(:,:,:,i-1),sst_extended(:,:,1:5,i));
end

% Insert last 5 days from previous year prior to the start of each year
for k = 1:30
    sst_long2(:,:,:,k) = cat(3,sst_extended(:,:,end-4:end,k),sst_long1(:,:,:,k));
end

% Calculate mean of each day within the 11-day window
% Set up array to store climatological mean
SST_Mean = [];
windowHalfWidth = 5;

for lon = 1:size(sst_long2,1)
    for lat = 1:size(sst_long2,2)
        for day = 6:371     % Equals to 1 Jan to 31 Dec
            SST_Mean(lon,lat,day-5) = mean(sst_long2(lon,lat,...
            day-windowHalfWidth:day+windowHalfWidth,:),'all','omitnan');
        end
    end
end

%% Calculate 90th percentile within the 11-day window
% Set up array to store climatological threshold
SST_pct = [];
for lon = 1:size(sst_long2,1)
    for lat = 1:size(sst_long2,2)
        for day = 6:371     % Equals to 1 Jan to 31 Dec
            SST_pct(lon,lat,day-5) = prctile(sst_long2(lon,lat,...
            day-windowHalfWidth:day+windowHalfWidth,:),90,'all');
        end
    end
end

%% Apply 31-day moving average 
sst_90long = smoothdata(cat(3,SST_pct,SST_pct,SST_pct),3,'movmean',31);

% Remove extra rows
SST_90 = sst_90long(:,:,367:367+365);

%% Select data from 1983-2012 only
SST_All = sst_extended(:,:,:,2:31);     

%% Calculate temperature anomaly
sst_anomaly = SST_All - SST_Mean;

% Separate leap and normal years
sst_anomaly_lp = sst_anomaly(:,:,:,2:4:end);
sst_anomaly_nm = sst_anomaly;
sst_anomaly_nm(:,:,:,2:4:end) = [];
sst_anomaly_nm(:,:,60,:) = [];  % Remove leap row

% Combined SST anomalies of all years into a 3D matrix
SST_AnomalyAll = [];

% Concatenate leap year in between normal years
k = [1:7; 2:3:22];
for i = 1:7
    sst_Anomaly1 = cat(3,sst_anomaly_lp(:,:,:,k(2*i-1)),sst_anomaly_nm(:,:,:,k(2*i)),...
    sst_anomaly_nm(:,:,:,k(2*i)+1),sst_anomaly_nm(:,:,:,k(2*i)+2));
    SST_AnomalyAll = cat(3,SST_AnomalyAll,sst_Anomaly1);
end

% Concatenate 1983 and 2012
SST_AnomalyAll = cat(3,sst_anomaly_nm(:,:,:,1),SST_AnomalyAll,sst_anomaly_lp(:,:,:,end));  

%% Detect SST above threshold, generate logical structure where 1 corresponds to SST above threshold
sst_warm = SST_All>SST_90;

% Separate leap and normal years
warm_lp = sst_warm(:,:,:,2:4:end);
warm_nm = sst_warm;
warm_nm(:,:,:,2:4:end) = [];
warm_nm(:,:,60,:) = [];     % Remove leap row

% Combined logical data of all years into a 3D matrix
SST_Warm = [];

% Concatenate leap year in between normal years
k = [1:7;2:3:22];
for i = 1:7
    sst_Warm1 = cat(3,warm_lp(:,:,:,k(2*i-1)),warm_nm(:,:,:,k(2*i)),...
    warm_nm(:,:,:,k(2*i)+1),warm_nm(:,:,:,k(2*i)+2));
    SST_Warm = cat(3,SST_Warm,sst_Warm1);
end

% Concatenate 1983 and 2012
SST_Warm = cat(3,warm_nm(:,:,:,1),SST_Warm,warm_lp(:,:,:,end));

%% Assign event numbers for adjacent SSTs above threshold
Z = false(nLon,nLat,1);             % Create a FALSE row
SST_Warm = cat(3,Z,SST_Warm);       % Insert false row as the first row to correct 'diff' error
SST_Warm2 = cumsum(cat(3,Z,(diff(SST_Warm,1,3))==1),3);
warm_evt = SST_Warm2.*SST_Warm;
warm_evt(~warm_evt) = nan;          % Change 0 to nan
warm_evt(:,:,1) = [];               % Remove the FALSE row, Z

%% Detect and count potential MHW event, when SSTs of one grid exceed
% threshold for >=5 consecutive days (for each cell, when the frequency of one event number >=5)
potential_evt = NaN(size(warm_evt));
for i = 1:nLon
    for j = 1:nLat
        n = 1;
        for count = 1:max(warm_evt(i,j,:),[],3,'omitnan')
            [~,~,D3] = findND(warm_evt(i,j,:)==count);            
            if max(D3)-min(D3)+1 >= 5
                potential_evt(i,j,min(D3):max(D3)) = n;
                n = n+1;
            end                
        end
    end
end

%% Clear un-needed variables 
clearvars -except SST_90 SST_All SST_Mean potential_evt nLon nLat SST_AnomalyAll YYYYmmdd YYYY lat_used lon_used 

%% Join events where gap<=2
MHW_Evt = potential_evt;

for i = 1:nLon
    for j = 1:nLat
        for count = 1:max(potential_evt(i,j,:),[],3,'omitnan')
            [~,~,A3] = findND(potential_evt(i,j,:)==count);
            [~,~,B3] = findND(potential_evt(i,j,:)==count+1);
            if ~isempty(A3) && ~isempty(B3)
                gap = min(B3)-max(A3)-1;        % Calculate duration of gap
                if gap <= 2
                    MHW_Evt(i,j,max(A3)+1:end) = MHW_Evt(i,j,max(A3)+1:end)-1;
                end                
            end                
        end
    end
end

%% Import intensity of each MHW from SST anomalies
MHW_Int = MHW_Evt;
for day = 1:size(MHW_Int,3)
    for lon = 1:nLon
        for lat = 1:nLat
            if MHW_Int(lon,lat,day) > 0         % if MHW exists
                MHW_Int(lon,lat,day) = SST_AnomalyAll(lon,lat,day);
            else 
                MHW_Int(lon,lat,day) = 0;
            end
        end
    end
end

%% Clear un-needed variables
clearvars -except nLon nLat YYYYmmdd YYYY lat_used lon_used MHW_Int MHW_Evt SST_90 SST_All SST_Mean MHW_All

%% Create array to record all MHW events
MHW_All=[];

for i = 1:nLon
    for j = 1:nLat
        for count = 1:max(MHW_Evt(i,j,:),[],3,'omitnan')
            if isfinite(count)
                [~,~,D3] = findND(MHW_Evt(i,j,:)==count);
                
                start_here = min(D3);
                end_here = max(D3);
                MHW_Start = YYYYmmdd(start_here);
                MHW_StartYr = YYYY(start_here);
                MHW_End = YYYYmmdd(end_here);
                MHW_EndYr = YYYY(end_here); 
                MHW_Duration = end_here-start_here+1;
                
                MHW_MaxInt = max(MHW_Int(i,j,D3),[],'omitnan');         % Maximum intensity
                MHW_CumInt = sum(MHW_Int(i,j,D3));                      % Cumulative intensity
                MHW_MeanInt = mean(MHW_Int(i,j,D3(1):D3(end)),'all');   % Mean intensity
                MHW_VarInt = std(MHW_Int(i,j,D3(1):D3(end)));           % Intensity variability
            
                x_loc = i;
                y_loc = j;
                longitude = lon_used(i);
                latitude = lat_used(j);
            
                MHW_All = [MHW_All;[MHW_Start MHW_End MHW_Duration MHW_MaxInt...
                    MHW_MeanInt MHW_VarInt MHW_CumInt MHW_StartYr MHW_EndYr longitude latitude]];
            end            
        end
    end
end 

%% Change array to table
MHW_table = [];
MHW_table = array2table(MHW_All, 'VariableNames', {'MHW_Start', 'MHW_End',...
    'Duration', 'MaxInt', 'MeanInt', 'VarInt', 'CumInt', 'MHW_StartYr', 'MHW_EndYr','Longitude','Latitude'});

%% Clear un-needed variables
clearvars -except SST_90 SST_All SST_Mean MHW_table MHW_All

%% Save MHW_All to directory path for further analysis
save("MHW_All","MHW_All");

%% Concatenate MHW_All if needed
% load MHW_All1 MHW_All2 MHW_All3   %% example filenames
% MHW_All = [MHW_All1;MHW_All2;MHW_All3];
