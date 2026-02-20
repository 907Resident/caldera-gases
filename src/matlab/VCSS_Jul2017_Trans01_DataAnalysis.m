
%% VCSS_Data Analysis

% VCSS_DataAnalysis
% Geoscientists: Ajayi, M. & Ayers, J.
% Location: Sulphur Springs, Valles Caldera, NM (Valles Caldera Sulphur
% Springs = VCSS)
% Data Type: Transect 01 | Pts. 01-05 (measured on 08 and 10 July 2017)
% Metadata File: (1) "...\Valles\08July2017\Picarro_Metadata_08Jul2017.txt"
%                (2) "...\Valles\10July2017\Picarro_Metadata_10Jul2017.txt"

%% Extra Functions Used in this Script
    % The following functions that were made by Ajayi, M. Their
    % documentation should be located in the "Functions" folder under
    % <"R:/Ayers/FrackingTN/Data/Matlab/Data Analysis">
    
        % lin_regress()
        % KeelingCurve()
        % StatsPlot_lite()
        % CH4_vs_Time()
        % CO2_vs_Time()
        % fullscreen()
        
%% Import Data Pre-Proccessed Data



% Date of measurements
ddmmmyyyy   = "07Jul2017";

% Data should be imported via "Merge_PicData_Fxn"
 directory  = working_dir+"*.dat";
[TT_PicData,PD_mtrx_nm,PD_mtrx_dt, dte] =                               ...
 Merge_PicData_Fxn(directory, 'America/Denver');

%% Gather the VCSS Data from the Forerunner Matrix 'PD_mtrx'

% Time Data (fraction of a day)
VCSS_X1         = VCSS_Time;                    % unitless
% Time Data (seconds relative to the start of the measuring period)
     % Change the serial time(fraction of a day) into seconds (the second
     % of the day, anywhere between 1 and 86400
VCSS_X1_sec     = 24 * 3600 .* VCSS_X1;     % seconds
    % Change the seconds to make them relative to t_0 (from 0 to 2400, or
    % 40 min)
VCSS_X1_sec     = VCSS_X1_sec - ...
                              min(VCSS_X1_sec);

% Important Metrics
    % These metrics have been organized in: import_VCSS_ALL_DATA

% Filter out Erroneous Measurements
VCSS_HP_12CH4_dry(VCSS_HP_12CH4_dry >= 1E06)            ...
                                                                     = NaN;
VCSS_HP_Delta_iCH4_Raw(VCSS_HP_Delta_iCH4_Raw >= 1E06)  ...
                                                                     = NaN;
VCSS_HR_12CH4_dry(VCSS_HR_12CH4_dry >= 1E06)            ...
                                                                     = NaN;
VCSS_Delta_iCO2_Raw(VCSS_Delta_iCO2_Raw >= 1E06)        ...
                                                                     = NaN;
VCSS_CO2_dry(VCSS_CO2_dry >= 1E06)                      ...
                                                                     = NaN;

%% Incorporate User Options

writeXL = input('Enter ''1'' to write main variables to an Excel file, enter zero if no Excel files are to be wrtitten: ');

%%  Get Rid of the Duplicate Measurement Values
    % In some instances, users of the Picarro have observed duplicate
    % measurements. That is, measurements over 1-3 seconds that are exactly
    % the same.  Thus, it becomes necessary to remove these duplicate
    % measurements.

% The removal of the duplicate measurements can be accomplished by creating
% a vector 1's and 0's.  The 1's will be included in rows that are kept and
% the 0's will indicate rows that are to be removed. 
    % For the VCSS, VCSS, and BHCS sites, the data presented groups of four
    % thus, only every fourth row is retainted
RN = zeros([length(VCSS_HP_Delta_iCH4_Raw), 1]);
    % Assigns '1' to every fourth row in the vector RN
for i = 1:4:length(RN)
RN(i) = 1;
end
    % Create a new matrix that will contain the desired measurements and
    % the keep/delete indicator by horizontally concatenating RN to the
    % measurement vectors
RN_CH4              = [RN, VCSS_X1, VCSS_HP_12CH4_dry, VCSS_HR_12CH4_dry...
                       VCSS_HP_Delta_iCH4_Raw, VCSS_HR_Delta_iCH4_Raw,  ...
                       VCSS_Delta_iCO2_Raw, VCSS_CO2_dry, VCSS_H2O,     ...
                       VCSS_Alarm];
    % Change the soon-to-be deleted rows to NaNs
for i = 1:length(RN)
    if RN_CH4(i,1) == 0
       RN_CH4(i,:) = NaN;
    end
    %*%* This line should be probably be deleted %*%*
    % VCSS_X1(i) = RN_CH4(i,2);
end
    % Delete the rows that contain NaNs (as determined by the indicator
    % vector)
RN_CH4(any(isnan(RN_CH4),2),:) = [];

% Filter out the extraneous date-timestamps
    idx         = 1:4:length(VCSS_Date);
    R_VCSS_Date = VCSS_Date(idx);

% Reassign the variables so that they represent the filtered information 
    VCSS_DateTime           = R_VCSS_Date;
    % Round all values to the beginning of the nearest second
        VCSS_DateTime = dateshift(VCSS_DateTime(:), 'start', 'second');
    VCSS_DateSerial         = datenum(VCSS_DateTime);
    VCSS_X1                 = RN_CH4(:,2);
    VCSS_X1_sec             = 24 * 3600 .* VCSS_X1;
    VCSS_X1_sec             = VCSS_X1_sec - ...
                                          min(VCSS_X1_sec);
    VCSS_HP_12CH4_dry       = RN_CH4(:,3);
    VCSS_HR_12CH4_dry       = RN_CH4(:,4);
    VCSS_HP_Delta_iCH4_Raw  = RN_CH4(:,5);
    VCSS_HR_Delta_iCH4_Raw  = RN_CH4(:,6);
    VCSS_Delta_Raw_iCO2_Raw = RN_CH4(:,7);
    VCSS_CO2_dry            = RN_CH4(:,8);
    VCSS_H2O                = RN_CH4(:,9);
    VCSS_Alarm              = RN_CH4(:,10);
    
%------Figure 01-------%
StatsPlot_lite(VCSS_DateSerial, VCSS_HP_12CH4_dry, 0)
    title('Sulphur Springs [CH_4] vs. Timestamp (ALL)')
    ylabel('[CH_4] (ppm)')
    legend('[CH_4]', '± 1 std', 'Mean', 'Max', 'Min', ...
           'Location', 'Northwest')
       
%------Figure 02-------%
StatsPlot_lite(VCSS_DateSerial, VCSS_CO2_dry, 13)
    title('Sulphur Springs [CO_2] vs. Timestamp (ALL)')
    ylabel('[CO_2] (ppm)')
    legend('[CO_2]', '± 1 std', 'Mean', 'Max', 'Min', ...
           'Location', 'Northwest')

clearvars R_VCSS_Date idx

%% Request the number of points on the transect

nchams = input('How many cycles of closed-open chamber measurements took place?: ');
            
%% Assign vectors for the time intervals for 'DateTime'

% Pre-background measurements
    % Pre-Background
        VCSS_PreBack_DT     = datetime(['08-Jul-2017 12:08:14';           
                                        '08-Jul-2017 12:50:00';           
                                        '08-Jul-2017 13:42:29';           
                                        '10-Jul-2017 10:57:30';           
                                        '10-Jul-2017 11:49:04';            
                                        '08-Jul-2017 12:14:33';           
                                        '08-Jul-2017 12:57:21';           
                                        '08-Jul-2017 13:47:53';           
                                        '10-Jul-2017 11:05:16';
                                        '10-Jul-2017 11:55:52']);
       VCSS_PreBack_DT      = reshape(VCSS_PreBack_DT, [nchams 2]);
       
    % Chamber ON
        VCSS_ChamON_DT      = datetime(['08-Jul-2017 12:14:37';
                                        '08-Jul-2017 12:57:24';
                                        '08-Jul-2017 13:47:57';
                                        '10-Jul-2017 11:05:20';
                                        '10-Jul-2017 11:55:45';
                                        '08-Jul-2017 12:27:09';
                                        '08-Jul-2017 13:03:04';
                                        '08-Jul-2017 14:06:55';
                                        '10-Jul-2017 11:24:54';
                                        '10-Jul-2017 12:03:02']);
        VCSS_ChamON_DT      = reshape(VCSS_ChamON_DT, [nchams 2]);
        
    % Post-Background
        VCSS_PostBack_DT    = datetime(['08-Jul-2017 12:30:20';
                                        '08-Jul-2017 13:04:33';
                                        '08-Jul-2017 14:08:38';
                                        '10-Jul-2017 11:25:50';
                                        '10-Jul-2017 12:03:06'; 
                                        '08-Jul-2017 12:35:21';
                                        '08-Jul-2017 13:09:46';
                                        '08-Jul-2017 14:11:35';
                                        '10-Jul-2017 11:30:39';
                                        '10-Jul-2017 12:03:31']);                                  
        VCSS_PostBack_DT      = reshape(VCSS_PostBack_DT, [nchams 2]);
                                  
    % Date Serials
        % Pre-Background
            VCSS_PreBack_DS     = datenum(VCSS_PreBack_DT);
        % Chamber ON
            VCSS_ChamON_DS      = datenum(VCSS_ChamON_DT);
        % Post-Background
            VCSS_PostBack_DS    = datenum(VCSS_PostBack_DT);

%%  Isolate [...]

% Pre-allocate vectors for pertinent metrics
    % Date
    VCSS_PreB_Date                      = NaT(350,  nchams);
    VCSS_ChamON_Date                    = NaT(350,  nchams);
    VCSS_PostB_Date                     = NaT(350,  nchams);
    % Seconds
    VCSS_PreB_sec                       = zeros(350, nchams);
    VCSS_ChamON_sec                     = zeros(350, nchams);
    VCSS_PostB_sec                      = zeros(350, nchams);
    % HP_12CH4_dry
    VCSS_PreB_HP_12CH4                  = NaN(350,  nchams);
    VCSS_ChamON_HP_12CH4                = NaN(350,  nchams);
    VCSS_PostB_HP_12CH4                 = NaN(350,  nchams);
    % HR_12CH4_dry
    VCSS_PreB_HR_12CH4                  = NaN(350, nchams);
    VCSS_ChamON_HR_12CH4                = NaN(350, nchams);
    VCSS_PostB_HR_12CH4                 = NaN(350, nchams);
    % HP_Delta_iCH4_Raw
    VCSS_PreB_HP_Delta_iCH4             = NaN(350, nchams);
    VCSS_ChamON_HP_Delta_iCH4           = NaN(350, nchams);
    VCSS_PostB_HP_Delta_iCH4            = NaN(350, nchams);
    % HR_Delta_iCH4_Raw
    VCSS_PreB_HR_Delta_iCH4             = NaN(350, nchams);
    VCSS_ChamON_HR_Delta_iCH4           = NaN(350, nchams);
    VCSS_PostB_HR_Delta_iCH4            = NaN(350, nchams);
    % CO2_dry
    VCSS_PreB_CO2                       = NaN(350, nchams);
    VCSS_ChamON_CO2                     = NaN(350, nchams);
    VCSS_PostB_CO2                      = NaN(350, nchams);
    % Delta_Raw_iCO2_Raw
    VCSS_PreB_Delta_iCO2                = NaN(350, nchams);
    VCSS_ChamON_Delta_iCO2              = NaN(350, nchams);
    VCSS_PostB_Delta_iCO2               = NaN(350, nchams);
    
for i = 1:nchams
% Gather the correct timestamps (serial numbers) to indicate the proper
% truncating locations in the dataset
PreB_counter_srt    = VCSS_PreBack_DS(i,1); 
PreB_counter_stp    = VCSS_PreBack_DS(i,2);

ChamON_counter_srt  = VCSS_ChamON_DS(i,1); 
ChamON_counter_stp  = VCSS_ChamON_DS(i,2);

PostB_counter_srt   = VCSS_PostBack_DS(i,1);
PostB_counter_stp   = VCSS_PostBack_DS(i,2);

% Use the find() function to acquire the proper indices to truncate the
% data
PreB_idx_srt    = find(PreB_counter_srt == VCSS_DateSerial);
PreB_idx_stp    = find(PreB_counter_stp == VCSS_DateSerial);

ChamON_idx_srt  = find(ChamON_counter_srt == VCSS_DateSerial);
ChamON_idx_stp  = find(ChamON_counter_stp == VCSS_DateSerial);

PostB_idx_srt   = find(PostB_counter_srt == VCSS_DateSerial);
PostB_idx_stp   = find(PostB_counter_stp == VCSS_DateSerial);
            
% Truncate the data
        %%% Date %%%
VCSS_PreB_Date(1:length(PreB_idx_srt:PreB_idx_stp),i)                   ...
                       = VCSS_DateTime(PreB_idx_srt:PreB_idx_stp);
VCSS_ChamON_Date(1:length(ChamON_idx_srt:ChamON_idx_stp),i)             ...
                       = VCSS_DateTime(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_Date(1:length(PostB_idx_srt:PostB_idx_stp),i)                ...
                       = VCSS_DateTime(PostB_idx_srt:PostB_idx_stp);
        %%% Time (Relative Sec) %%%
    % See 'Seconds' below                           
        %%% HP [12CH4] %%%
VCSS_PreB_HP_12CH4(1:length(PreB_idx_srt:PreB_idx_stp),i)               ...
                       =  VCSS_HP_12CH4_dry(PreB_idx_srt:PreB_idx_stp);
VCSS_ChamON_HP_12CH4(1:length(ChamON_idx_srt:ChamON_idx_stp),i)         ...
                       =  VCSS_HP_12CH4_dry(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_HP_12CH4(1:length(PostB_idx_srt:PostB_idx_stp),i)            ...
                       =  VCSS_HP_12CH4_dry(PostB_idx_srt:PostB_idx_stp);
       %%% HR [12CH4] %%%
VCSS_PreB_HR_12CH4(1:length(PreB_idx_srt:PreB_idx_stp),i)               ...
                       =  VCSS_HR_12CH4_dry(PreB_idx_srt:PreB_idx_stp);
VCSS_ChamON_HR_12CH4(1:length(ChamON_idx_srt:ChamON_idx_stp),i)         ...
                       =  VCSS_HR_12CH4_dry(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_HR_12CH4(1:length(PostB_idx_srt:PostB_idx_stp),i)            ...
                       =  VCSS_HR_12CH4_dry(PostB_idx_srt:PostB_idx_stp);
     %%% HP d13C-CH4 %%% 
VCSS_PreB_HP_Delta_iCH4(1:length(PreB_idx_srt:PreB_idx_stp),i)          ...
                       =  VCSS_HP_Delta_iCH4_Raw(PreB_idx_srt:PreB_idx_stp);     
VCSS_ChamON_HP_Delta_iCH4(1:length(ChamON_idx_srt:ChamON_idx_stp),i)    ...
                       =  VCSS_HP_Delta_iCH4_Raw(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_HP_Delta_iCH4(1:length(PostB_idx_srt:PostB_idx_stp),i)       ...
                       =  VCSS_HP_Delta_iCH4_Raw(PostB_idx_srt:PostB_idx_stp);
     %%% HR d13C-CH4 %%%
VCSS_PreB_HR_Delta_iCH4(1:length(PreB_idx_srt:PreB_idx_stp),i)          ...
                       =  VCSS_HR_Delta_iCH4_Raw(PreB_idx_srt:PreB_idx_stp);     
VCSS_ChamON_HR_Delta_iCH4(1:length(ChamON_idx_srt:ChamON_idx_stp),i)    ...
                       =  VCSS_HR_Delta_iCH4_Raw(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_HR_Delta_iCH4(1:length(PostB_idx_srt:PostB_idx_stp),i)       ...
                       =  VCSS_HR_Delta_iCH4_Raw(PostB_idx_srt:PostB_idx_stp);
    %%% [CO2_dry] %%%
VCSS_PreB_CO2(1:length(PreB_idx_srt:PreB_idx_stp),i)                    ...
                       =  VCSS_CO2_dry(PreB_idx_srt:PreB_idx_stp);     
VCSS_ChamON_CO2(1:length(ChamON_idx_srt:ChamON_idx_stp),i)              ...
                       =  VCSS_CO2_dry(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_CO2(1:length(PostB_idx_srt:PostB_idx_stp),i)                 ...
                       =  VCSS_CO2_dry(PostB_idx_srt:PostB_idx_stp);
    %%% Delta_Raw_iCO2_Raw %%%
VCSS_PreB_Delta_iCO2(1:length(PreB_idx_srt:PreB_idx_stp),i)             ...
                       =  VCSS_Delta_Raw_iCO2_Raw(PreB_idx_srt:PreB_idx_stp);    
VCSS_ChamON_Delta_iCO2(1:length(ChamON_idx_srt:ChamON_idx_stp),i)       ...
                       =  VCSS_Delta_Raw_iCO2_Raw(ChamON_idx_srt:ChamON_idx_stp);
VCSS_PostB_Delta_iCO2(1:length(PostB_idx_srt:PostB_idx_stp),i)          ...
                       =  VCSS_Delta_Raw_iCO2_Raw(PostB_idx_srt:PostB_idx_stp);                   
end

            
% Seconds
    %PreB
        j = 1;
    while j <= 5
        for i = 1:length(VCSS_ChamON_Date)-1
            % Pre-Background
            s_PreB                  = diff(VCSS_PreB_Date(:,j));
            ds_PreB                 = seconds(s_PreB);
            VCSS_PreB_sec(i+1,j)  = VCSS_PreB_sec(i,j) + ds_PreB(i);
        end
        j                           = j + 1;
    end
    
    % ChamON
    j = 1;
    while j <= 5
        for i = 1:length(VCSS_ChamON_Date)-1
            % Chamber ON
            s_ChamON                = diff(VCSS_ChamON_Date(:,j));
            ds_ChamON               = seconds(s_ChamON);
            VCSS_ChamON_sec(i+1,j)  = VCSS_ChamON_sec(i,j) + ds_ChamON(i);
        end
        j                           = j + 1;
    end
    
    %PostB
        j = 1;
    while j <= 5
        for i = 1:length(VCSS_ChamON_Date)-1
            % Post-Background
            s_PostB                 = diff(VCSS_PostB_Date(:,j));
            ds_PostB                 = seconds(s_PostB);
            VCSS_PostB_sec(i+1,j)  = VCSS_PreB_sec(i,j) + ds_PreB(i);
        end
        j                           = j + 1;
    end
    
    % Return extraneous zeros to NaN (this will help in subsequent analyses
    VCSS_PreB_sec(VCSS_PreB_sec(2:end) == 0)                         = NaN;
    VCSS_ChamON_sec(VCSS_ChamON_sec(2:end) == 0)                     = NaN;
    VCSS_PostB_sec(VCSS_PostB_sec(2:end) == 0)                       = NaN;
    
%% Preliminary Visualization of the ISMs    
    
% Concatenate all of the ISMs into a three dimensional array, where
% the third dimension represents the individual chamber measurements
VCSS_ISM_DATA = zeros([length(VCSS_ChamON_HP_12CH4), ...
                                    7, nchams]);

% NOTE: at this point there are only six ISM loaded into the matrix at this
% point. The remaining ISMs (the concentrations in ug/L) will be introduced
% at a later point in the script.
for idx = 1:nchams
        VCSS_ISM_DATA(:,1,idx) = VCSS_ChamON_sec(:,idx);
        VCSS_ISM_DATA(:,1,idx) = VCSS_ISM_DATA(:,1,idx) ./ hour;
        VCSS_ISM_DATA(:,2,idx) = VCSS_ChamON_HP_12CH4(:,idx);
        VCSS_ISM_DATA(:,3,idx) = VCSS_ChamON_HR_12CH4(:,idx);
        VCSS_ISM_DATA(:,4,idx) = VCSS_ChamON_HP_Delta_iCH4(:,idx);
        VCSS_ISM_DATA(:,5,idx) = VCSS_ChamON_HR_Delta_iCH4(:,idx);
        VCSS_ISM_DATA(:,6,idx) = VCSS_ChamON_CO2(:,idx);
        VCSS_ISM_DATA(:,7,idx) = VCSS_ChamON_Delta_iCO2(:,idx);
end

% Preliminary ISM Data Visualizations
% Concentration-Time Plots
%-----Figure 04-----%
% Create a 8 x nchams/2 array scatter plots of concentration vs time (rel
% sec) for both gases
    fullscreen,
for idx = 1:nchams
    
%--CH4--%
    % This if-else block helps to determine where to place each of the
    % graphs 
    if idx <= 3 
        CH4_pl = idx;
    else
        CH4_pl = idx+3;
    end
    subplot(4,ceil(nchams/2),CH4_pl)
    s_CH4 = scatter(VCSS_ChamON_sec(:,idx),VCSS_ISM_DATA(:,2,idx));
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
            if idx     <= 5
                trans   = 1;
                pnt     = idx;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       'VCSS %d.%d Jul2017 [CH_4] vs. Rel Time',      ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
            % The maximum value on the x-axis will be 1200 (= 60 sec *
            % 20 min)
            xlim([0 1200])
            % Provide tick marks in 300 sec (5 min) intervals
            xticks(0:300:1200)
            xlabel('Relative Time (sec)')
            % Apply gridding on both axis to see the data better
            grid on
            
%--CO2--%
    % This if-else block helps to determine where to place each of the
    % graphs
    if idx <= 3 
        CO2_pl = idx+3;
    else
        CO2_pl = idx+6;
    end 
    subplot(4,ceil(nchams/2),CO2_pl)
    s_CO2 = scatter(VCSS_ChamON_sec(:,idx),VCSS_ISM_DATA(:,6,idx));
    s_CO2.MarkerEdgeColor = 'm';
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
            if idx     <= 5
                trans   = 1;
                pnt     = idx;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       'VCSS %d.%d Jul2017 [CO_2] vs. Rel Time',      ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
            % The maximum value on the x-axis will be 1200 (= 60 sec *
            % 20 min)
            xlim([0 1200])
            % Provide tick marks in 300 sec (5 min) intervals
            xticks(0:300:1200)
            xlabel('Relative Time (sec)')
            % Apply gridding on both axis to see the data better
            grid on        
                
end

%% Convert gas measurements from ppm tg cgs units - **Unexcavated**

% Converting constants
Vm              = 22.71108; % Standard Molar Volume (L)
CH4_molec_wt    = 16.04;    % g/mol
CO2_molec_wt    = 44.01;    % g/mol

% Imported concentrations (resampled)
c_CH4_ppm_UnEx  = VCSS_ISM_DATA(:,2,:);   
c_CO2_ppm_UnEx  = VCSS_ISM_DATA(:,6,:);
% Imported times (resampled)
t_UnEx          = VCSS_ISM_DATA(:,1,:);
% Height (m)
    % 9.48 cm = 0.0948 m
H = 0.0948;

%%% To get fluxes, think of concentration expressed in milligrams per cubic
%%% meter mg/m3
    %%% Convert from ratio of gas to air (ppm) to mass per volume (mg /
    %%% m^3)
        %%% (Molec_wt / Vm) * ppm = mg / m^3
        c_CH4_cgs_UnEx = (CH4_molec_wt / Vm) .* c_CH4_ppm_UnEx;    % mg m-3
    %%% Convert from mg / L to mg / m3
        %%% 1 ug/L = 1000 ug/m3 %%% 
%         c_CH4_cgs_UnEx =  1000 .* c_CH4_cgs_UnEx;
    %%% Repeat for carbon dioxide     
        c_CO2_cgs_UnEx = (CO2_molec_wt / Vm) .* c_CO2_ppm_UnEx;    % mg m-3
%         c_CO2_cgs_UnEx =  1000 .* c_CO2_cgs_UnEx;
        
 for idx = 1:nchams
        VCSS_ISM_DATA(:,8,idx) = c_CH4_cgs_UnEx(:,1,idx);
        VCSS_ISM_DATA(:,9,idx) = c_CO2_cgs_UnEx(:,1,idx);
end       
        
%% Flux Calculation

%----||Linear Model Analysis||----%

    % A two-minute equilibration period (approx 33 measurements) is allowed
    % before an assessment of the linear fit.  After this two-minute
    % period, the next three minutes are analyzed (approx. 52 measurements)
    % using a simple linear regression model
        % C_t = b + m*t, where the slope, m, is dC_dt
        
        %----Methane----%
    % Pre-allocate vector for fluxes (per measurement cycle)
    CH4_lin_mdl     = cell ([1 nchams]); 
    CH4_lin_slope   = zeros([2 4 nchams]);
    CH4_lin_flux    = zeros([1 nchams]);
    % Use for-loop to run a linear regression for each of the chamber
    % measuerments
    for i = 1:nchams
       CH4_lin_mdl{i}           = fitlm(VCSS_ISM_DATA(:,1,i),           ...
                                        VCSS_ISM_DATA(:,8,i));
       CH4_lin_slope(:,:,i)     = table2array(CH4_lin_mdl{i}.Coefficients);
       % Quantify the flux for each measurement (mg m^-2 hr^-1)
            % The time component has already been converted to hours
       CH4_lin_flux(i)          = CH4_lin_slope(2,1,i) .* H;
    end
    
        %----Carbon Dioxide----%
    % Pre-allocate vector for fluxes (per measurement cycle)
    CO2_lin_mdl     = cell ([1 nchams]); 
    CO2_lin_slope   = zeros([2 4 nchams]);
    CO2_lin_flux    = zeros([1 nchams]);
    % Use for-loop to run a linear regression for each of the chamber
    % measuerments
    for i = 1:nchams
       CO2_lin_mdl{i}           = fitlm(VCSS_ISM_DATA(:,1,i), ...
                                        VCSS_ISM_DATA(:,9,i));
       CO2_lin_slope(:,:,i)     = table2array(CO2_lin_mdl{i}.Coefficients);
       % Quantify the flux for each measurement (mg m^-2 hr^-1)
       CO2_lin_flux(i)          = CO2_lin_slope(2,1,i) .* H;
    end         
        
%----||Exponential Model Analysis||----%
    
    %----Methane----%

    % A two-minute equilibration period (approx 33 measurements) is allowed
    % before an assessment of the exponential fit.  After this two-minute
    % period, the remaining portion of the measurement period is evaluted
    % using an exponential fit:
        % C_t = psi + (C_0 - psi) * exp(-kappa*t)
    % Pre-allocate vector for fluxes (per measurement cycle)
    CH4_exp_mdl     = cell([1 nchams]);
    % Define the options for the exponential fit (i.e. fitting method, starting
    % point for the iteration)
    CH4_exp_ft_ops  = fitoptions('StartPoint', [0 0], 'Method', 'NonlinearLeastSquares');
    % Create an exponential model to evaluate the flux data using fittype()
    CH4_exp_fxn     = fittype(@(psi, kappa, c_0, t) psi + (c_0 - psi) * exp(-kappa*t),  ...
                          'dependent', {'C_t'}, 'independent', {'t'},               ...
                          'coefficients', {'psi', 'kappa'}, 'problem', 'c_0',       ...
                          'options', CH4_exp_ft_ops);
    %c_0         = 1.3575; % mg / m^3 
        % (This is equivalent to 1.9 ppm CH4 and could be adjusted to
        % indicate the concentration at time = 0)
    CH4_exp_slope   = zeros([1 3 nchams]);
    CH4_exp_flux    = zeros([1 nchams]);

    x_expData = cell([1 nchams]);
    y_expData = cell([1 nchams]);

for i = 1:nchams
% Apply an exponential model to evaluate the flux data using fit()
    % The row indicies are hard coded to ensure that no NaN values are
    % analyzed; additionally, this segment of data accounts for the
    % majority of the data after the equilibriation period
[x_exp_prepData, y_exp_prepData]  = ...
                        prepareCurveData(VCSS_ISM_DATA(:,1,i), ...
                                         VCSS_ISM_DATA(:,8,i));   
x_expData{i}            = x_exp_prepData;
y_expData{i}            = y_exp_prepData;
CH4_exp_mdl{i}          = fit(x_expData{i},y_expData{i},CH4_exp_fxn,              ...
                             'problem', VCSS_ISM_DATA(1,8,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
CH4_exp_slope(1,1,i)  = CH4_exp_mdl{i}.psi;
CH4_exp_slope(1,2,i)  = CH4_exp_mdl{i}.kappa;
CH4_exp_slope(1,3,i)  = CH4_exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
CH4_exp_flux(i)       =((CH4_exp_mdl{i}.psi - CH4_exp_mdl{i}.c_0) .* CH4_exp_mdl{i}.kappa) .* H;
end 

clearvars x_expData y_expData x_exp_prepData y_exp_prepData

    %----Carbon Dioxide----%

    % A two-minute equilibration period (approx 33 measurements) is allowed
    % before an assessment of the exponential fit.  After this two-minute
    % period, the remaining portion of the measurement period is evaluted
    % using an exponential fit:
        % C_t = psi + (C_0 - psi) * exp(-kappa*t)
    % Pre-allocate vector for fluxes (per measurement cycle)
    CO2_exp_mdl     = cell([1 nchams]);
    % Define the options for the exponential fit (i.e. fitting method, starting
    % point for the iteration)
    CO2_exp_ft_ops  = fitoptions('StartPoint', [0 0], 'Method', 'NonlinearLeastSquares');
    % Create an exponential model to evaluate the flux data using fittype()
    CO2_exp_fxn     = fittype(@(psi, kappa, c_0, t) psi + (c_0 - psi) * exp(-kappa*t),  ...
                          'dependent', {'C_t'}, 'independent', {'t'},               ...
                          'coefficients', {'psi', 'kappa'}, 'problem', 'c_0',       ...
                          'options', CO2_exp_ft_ops);
    %c_0         = 1.3575; % mg / m^3 
        % (This is equivalent to 1.9 ppm CH4 and could be adjusted to
        % indicate the concentration at time = 0)
    CO2_exp_slope   = zeros([1 3 nchams]);
    CO2_exp_flux    = zeros([1 nchams]);

    x_expData = cell([1 nchams]);
    y_expData = cell([1 nchams]);

for i = 1:nchams
% Apply an exponential model to evaluate the flux data using fit()
    % The row indicies are hard coded to ensure that no NaN values are
    % analyzed; additionally, this segment of data accounts for the
    % majority of the data after the equilibriation period
[x_exp_prepData, y_exp_prepData]  = ...
                        prepareCurveData(VCSS_ISM_DATA(:,1,i), ...
                                         VCSS_ISM_DATA(:,9,i));   
x_expData{i}            = x_exp_prepData;
y_expData{i}            = y_exp_prepData;
CO2_exp_mdl{i}          = fit(x_expData{i},y_expData{i},CO2_exp_fxn,              ...
                             'problem', VCSS_ISM_DATA(1,9,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
CO2_exp_slope(1,1,i)  = CO2_exp_mdl{i}.psi;
CO2_exp_slope(1,2,i)  = CO2_exp_mdl{i}.kappa;
CO2_exp_slope(1,3,i)  = CO2_exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
CO2_exp_flux(i)       =((CO2_exp_mdl{i}.psi - CO2_exp_mdl{i}.c_0) .* CO2_exp_mdl{i}.kappa) .* H;
end 
clearvars x_expData y_expData x_exp_prepData y_exp_prepData

%% Isotope Data Anlysis

%----Methane----%

% Pre-allocate vectors for KeelingCurve Results
CH4_KC_rslts    = cell([1 nchams]);
CH4_KC_coeffs   = zeros([1 2 nchams]);
CH4_KC_coeffs_ci= zeros([2 2 nchams]);

for i = 1:nchams
% Prepare data for the regression (i.e. remove NaN values, shape the
% vectors as vertical vectors, etc.) 
[xData, yData] = prepareCurveData((1 ./ VCSS_ISM_DATA(:,2,i)), VCSS_ISM_DATA(:,4,i));

% Set up fittype and options.
KC_ft = fittype('poly1');

% Fit model to data.
[KC_fitresult, KC_gof] = fit(xData, yData, KC_ft);

% Gather the coefficients for the fit equation
fitvals = coeffvalues(KC_fitresult);
% Gather the confidence intervals for c
    CH4_KC_rslts{i}             = KC_fitresult;
    CH4_KC_coeffs(1,1,i)        = fitvals(1);
    CH4_KC_coeffs(1,2,i)        = fitvals(2);
    CH4_KC_coeffs_ci(:,:,i)     = confint(KC_fitresult, 0.95);
end
clearvars fitvals KC_ft xData yData KC_fitresult

%----Carbon Dioxide----%
CO2_KC_rslts    = cell([1 nchams]);
CO2_KC_coeffs   = zeros([1 2 nchams]);
CO2_KC_coeffs_ci= zeros([2 2 nchams]);

for i = 1:nchams
% Prepare data for the regression (i.e. remove NaN values, shape the
% vectors as vertical vectors, etc.) 
[xData, yData] = prepareCurveData((1 ./ VCSS_ISM_DATA(:,2,i)), VCSS_ISM_DATA(:,7,i));

% Set up fittype and options.
KC_ft = fittype('poly1');

% Fit model to data.
[KC_fitresult, KC_gof] = fit(xData, yData, KC_ft);

% Gather the coefficients for the fit equation
fitvals = coeffvalues(KC_fitresult);
% Gather the confidence intervals for c
    CO2_KC_rslts{i}             = KC_fitresult;
    CO2_KC_coeffs(1,1,i)        = fitvals(1);
    CO2_KC_coeffs(1,2,i)        = fitvals(2);
    CO2_KC_coeffs_ci(:,:,i)     = confint(KC_fitresult, 0.95);
end
clearvars fitvals KC_ft xData yData KC_fitresult

% Calculate alpha
    % ? = [(?13C-CO2 - ?13C-CH4)/1000 ] + 1
    alpha           = zeros([1,nchams]);
    thous_ln_alpha  = zeros([1,nchams]);
    for i = 1:nchams
    alpha(i)            = ((CO2_KC_coeffs(1,2,i) - CH4_KC_coeffs(1,2,i)) ./ 1E03) + 1;
    thous_ln_alpha(i)   = 1E03 .* log(alpha(i));
    end

%% Data Visualization

%----Methane----%
% Plot regressions with 95% confidence intervals for the slope   
fullscreen,
minimum = min(min(VCSS_ISM_DATA(:,8,:)));
maximum = max(max(VCSS_ISM_DATA(:,8,:)));
for i = 1:nchams
[y, yci]            = predict(CH4_lin_mdl{i}, VCSS_ISM_DATA(:,1,i), 'Simultaneous', 'on');
coeff_slope         = table2array(CH4_lin_mdl{1,i}.Coefficients(2,1));
coeff_incp          = table2array(CH4_lin_mdl{1,i}.Coefficients(1,1));
r_squared           = CH4_lin_mdl{1,i}.Rsquared.Ordinary;
coeff_ci            = coefCI(CH4_lin_mdl{i});
regress_ttl_string  = sprintf('VCSS | Lin. Model | Trans 01 | Pt. %d', i);
fit_output          = sprintf(' [CH_4] = %2.3E * t + %2.3E \n Slope CI (%2.3f, %2.3f) \n Y-Intercept CI (%2.3f, %2.3f) \n R-squared = %1.4f', ...
                      coeff_slope, coeff_incp,                          ...
                      coeff_ci(1), coeff_ci(2),                         ...
                      coeff_ci(3), coeff_ci(4),                         ...
                      r_squared);
subplot(3,2,i)    
    plot(VCSS_ISM_DATA(:,1,i), y, 'r');
        hold on, 
    plot(VCSS_ISM_DATA(:,1,i), yci(:,1), '--r'); 
        hold on, 
    plot(VCSS_ISM_DATA(:,1,i), yci(:,2), '--r'); 
        hold on, 
    plot(VCSS_ISM_DATA(:,1,i), VCSS_ISM_DATA(:,8,i), 'o'); 
        grid on
        title(regress_ttl_string, 'FontSize', 8)
        xlabel('Time (hr)', 'FontSize', 7)
        ylabel('[CH_4] (mg)', 'FontSize', 7)
        ylim([minimum maximum])
        xlim([0 0.25])
    if      i == 1
       annotation('textbox', [0.308 0.714 0.150 0.111], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 2
       annotation('textbox', [0.754 0.709 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 3
        annotation('textbox', [0.127 0.516 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif i == 4
        annotation('textbox', [0.568 0.518 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif i == 5
        annotation('textbox', [0.324 0.117 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    end
end

clearvars coeff_slope coeff_incp r_squared coeff_ci regress_ttl_string fit_output maximum minimum

%----Carbon Dioxide----%
% Plot regressions with 95% confidence intervals for the slope   
fullscreen,
minimum = min(min(VCSS_ISM_DATA(:,9,:)));
maximum = max(max(VCSS_ISM_DATA(:,9,:)));
for i = 1:nchams
[y, yci]            = predict(CO2_lin_mdl{i}, VCSS_ISM_DATA(:,1,i), 'Simultaneous', 'on');
coeff_slope         = table2array(CO2_lin_mdl{1,i}.Coefficients(2,1));
coeff_incp          = table2array(CO2_lin_mdl{1,i}.Coefficients(1,1));
r_squared           = CO2_lin_mdl{1,i}.Rsquared.Ordinary;
coeff_ci            = coefCI(CO2_lin_mdl{i});
regress_ttl_string  = sprintf('VCSS | Lin. Model | Trans 01 | Pt. %d', i);
fit_output          = sprintf(' [CO_2] = %2.3E * t + %2.3E \n Slope CI (%2.3f, %2.3f) \n Y-Intercept CI (%2.3f, %2.3f) \n R-squared = %1.4f', ...
                      coeff_slope, coeff_incp,                          ...
                      coeff_ci(1), coeff_ci(2),                         ...
                      coeff_ci(3), coeff_ci(4),                         ...
                      r_squared);
subplot(3,2,i)    
    plot(VCSS_ISM_DATA(:,1,i), y, 'r');
        hold on, 
    plot(VCSS_ISM_DATA(:,1,i), yci(:,1), '--r'); 
        hold on, 
    plot(VCSS_ISM_DATA(:,1,i), yci(:,2), '--r'); 
        hold on, 
    plot(VCSS_ISM_DATA(:,1,i), VCSS_ISM_DATA(:,9,i), 'mo'); 
        grid on
        title(regress_ttl_string, 'FontSize', 8)
        xlabel('Time (hr)', 'FontSize', 7)
        ylabel('[CH_4] (ppm)', 'FontSize', 7)
        ylim([minimum maximum])
        xlim([0 0.25])
    if      i == 1
       annotation('textbox', [0.308 0.714 0.150 0.111], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 2
       annotation('textbox', [0.754 0.709 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 3
        annotation('textbox', [0.127 0.516 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif i == 4
        annotation('textbox', [0.568 0.518 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif i == 5
        annotation('textbox', [0.324 0.117 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    end
end
clearvars coeff_slope coeff_incp r_squared coeff_ci regress_ttl_string fit_output maximum minimum y yci
    
    %-----Exponential Model-----%
%     e_flux_time = cell([length(CH4_exp_flux), 1]);
%     T_exp_flux          = table;
%     for i = 1:npts
%         flujo_tiempo    = datestr(VCSS_ISM_DATA(490,1,i) + datenum(2017,03,06), 0);
%         e_flux_time{i}  = datetime(flujo_tiempo,'InputFormat','dd-MM-yyyy HH:mm:ss');
%     end
% 
%     T_exp_flux.Time     = [e_flux_time{:,1}]';
%     T_exp_flux.Flux     = CH4_exp_flux(1,:)';
%     
% fullscreen,
%     plot(T_exp_flux.Time, T_exp_flux.Flux, 'LineStyle', '--',           ...
%         'Marker', 'd');
%             title('CH_{4} Fluxes (Lin. Model) | BHCS')
%             ylabel('Flux (mg m^{-2} hr^{-1}')
%     
% % Intercepts vs. Time
% 
%     incp = zeros([npts, 1]);
% for i = 1:npts
%    incp(i) = CH4_KC_coeffs(1,2,i); 
% end
% 
% T_incp = table;
%     T_incp.Time = [e_flux_time{:,1}]';
%     T_incp.Incp = incp;
%         clearvars incp
% 
% fullscreen,
%     plot(T_incp.Time, T_incp.Incp, 'LineStyle', '--', 'Marker', 'X', ...
%          'MarkerSize', 14);
%             title('Isotopic Composition of Methane | BHCS','Interpreter', 'tex')
%             ylabel('\delta^{13}C-CH_{4} (‰)','Interpreter', 'tex')
%             % Add the error at each point after the fact with textboxes           
            
% CO2/CH4 Ratio per location
CO2_CH4_ratio = zeros([length(VCSS_ISM_DATA), 1]);
for i = 1:nchams
    CO2_CH4_ratio(:,i) =  VCSS_ISM_DATA(:,6,i)./VCSS_ISM_DATA(:,2,i);
end

% CO2/CH4 Ratio v. Time
fullscreen,
minimum = min(min(CO2_CH4_ratio));
maximum = max(max(CO2_CH4_ratio));
for i = 1:nchams
    title_string = sprintf('VCSS | CO_2:CH_4 Ratio vs. Time | Trans 01 | Pt %d', i);
    subplot(3,2,i)
        scatter(VCSS_ISM_DATA(:,1,i), CO2_CH4_ratio(:,i));
        title(title_string)
        ylabel('[CO_2]:[CH_4] Ratio')
        xlabel('Time (hr)')
        ylim([minimum maximum]);
        xlim([0 0.25])
end
clearvars title_string maximum minimum

% CH4:CO2 Flux Ratio (mg m-2 hr-1)
    % Linear Model
Flx_Ratio = CH4_lin_flux ./ CO2_lin_flux;
    % Visualize
    fullscreen,
    subplot(1,2,1)
    b1 = bar(Flx_Ratio, 'FaceColor', [0.3010 0.7450 0.9330],            ...
            'EdgeColor', 'k');
            title('CH_4-CO_2 Flux Ratio | VCSS | Trans 01')
            ylabel('F_{CH_4}:F_{CO_2}')
            xticklabels({'Point 01', 'Point 02', 'Point 03', 'Point 04',...
                         'Point 05'});
    subplot(1,2,2)
    b2 = bar(log10(Flx_Ratio), 'FaceColor', [0.3010 0.7450 0.9330],            ...
            'EdgeColor', 'k');
            ylabel('log_{10}(F_{CO2_2}:F_{CH_4}) (mg m^{-2} hr^{-1})')
            xticklabels({'Point 01', 'Point 02', 'Point 03', 'Point 04',...
                         'Point 05'});
                     
% Compare carbon isotopes for both gases
VCSS_isocomp    = zeros([1 nchams]);
for i = 1:nchams
    VCSS_isocomp(i) = isocomp(VCSS_ISM_DATA(:,4,i),VCSS_ISM_DATA(:,7,i));
end

% Meteorological Data
    % Import Data
Import_MeteoData_08Jul2017
    % Visualize Meteorological Data
fullscreen,
subplot(3,1,1)
    bar(MeteoMeas_08Jul2017.MeteoTime, MeteoMeas_08Jul2017.AirTemp)
        title('Ambient Temperature')
        ylabel('Temperature (°C)')
subplot(3,1,2)
    bar(MeteoMeas_08Jul2017.MeteoTime, MeteoMeas_08Jul2017.BaroPress)
        title('Barometric Pressure')
        ylabel('Pressure (mbar)')
subplot(3,1,3)
    bar(MeteoMeas_08Jul2017.MeteoTime, MeteoMeas_08Jul2017.Precip)
        title('Total Precipitation')
        ylabel('Precipitation (mm)')

% Soil Measurements
    % TO BE COMPLETED [15 Feb 2018]

%% Export Data
    
% Write the chamber data into an Excel file
if writeXL == 1
    for i = 1:nchams
        export_chamON_data  = [VCSS_ISM_DATA(:,1,i), VCSS_ISM_DATA(:,2,i),   ...
                               VCSS_ISM_DATA(:,3,i), VCSS_ISM_DATA(:,4,i),   ...
                               VCSS_ISM_DATA(:,5,i), VCSS_ISM_DATA(:,6,i),   ...
                               VCSS_ISM_DATA(:,7,i), VCSS_ISM_DATA(:,8,i),   ...
                               VCSS_ISM_DATA(:,9,i)];
        export_tabname          = sprintf('VCSS_Trans01_Pt%d', i);
        export_lin_mdl_CH4      = CH4_lin_slope(:,:,i);
        export_lin_flux_CH4     = CH4_lin_flux(i);
        export_exp_mdl_CH4      = CH4_exp_slope(:,:,i);
        export_exp_flux_CH4     = CH4_exp_flux(i);
        export_KC_CH4_coeffs    = CH4_KC_coeffs(:,:,i);
        export_KC_CH4_coeffs_ci = CH4_KC_coeffs_ci(:,:,i);
        export_lin_mdl_CO2      = CO2_lin_slope(:,:,i);
        export_lin_flux_CO2     = CO2_lin_flux(i);
        export_exp_mdl_CO2      = CO2_exp_slope(:,:,i);
        export_exp_flux_CO2     = CO2_exp_flux(i);
        export_KC_CO2_coeffs    = CO2_KC_coeffs(:,:,i);
        export_KC_CO2_coeffs_ci = CO2_KC_coeffs_ci(:,:,i);
        export_alpha            = alpha(i);
        export_thous_ln_alpha   = thous_ln_alpha(i);
        XL_filename             = 'C:/Users/moyoa/Google Drive/CompSci/PhD Dissertation/Data Analysis/Picarro/Valles/VCSS_Trans01_Jul2017.xlsx';
            % Currently, when using the xlswrite function in a for-loop, a
            % new (blank) workbook is created on every iteration even
            % though the data is properly written into BHCS_Mar2017
        XL_Cham_ON          = xlswrite(XL_filename, export_chamON_data,        export_tabname, 'A3' );
        XL_lin_mdl_CH4      = xlswrite(XL_filename, export_lin_mdl_CH4,        export_tabname, 'K3' );
        XL_lin_flux_CH4     = xlswrite(XL_filename, export_lin_flux_CH4,       export_tabname, 'P3' );
        XL_exp_mdl_CH4      = xlswrite(XL_filename, export_exp_mdl_CH4,        export_tabname, 'R3' );
        XL_exp_flux_CH4     = xlswrite(XL_filename, export_exp_flux_CH4,       export_tabname, 'V3' );
        XL_KC_CH4_coeffs    = xlswrite(XL_filename, export_KC_CH4_coeffs,      export_tabname, 'X3' );
        XL_KC_CH4_coeffs_ci = xlswrite(XL_filename, export_KC_CH4_coeffs_ci,   export_tabname, 'Z3' );
        XL_lin_mdl_CO2      = xlswrite(XL_filename, export_lin_mdl_CO2,        export_tabname, 'AC3');
        XL_lin_flux_CO2     = xlswrite(XL_filename, export_lin_flux_CO2,       export_tabname, 'AH3');
        XL_exp_mdl_CO2      = xlswrite(XL_filename, export_exp_mdl_CO2,        export_tabname, 'AJ3');
        XL_exp_flux_CO2     = xlswrite(XL_filename, export_exp_flux_CO2,       export_tabname, 'AN3');
        XL_KC_CO2_coeffs    = xlswrite(XL_filename, export_KC_CO2_coeffs,      export_tabname, 'AP3');
        XL_KC_CO2_coeffs_ci = xlswrite(XL_filename, export_KC_CO2_coeffs_ci,   export_tabname, 'AR3');
        XL_alpha            = xlswrite(XL_filename, export_alpha,              export_tabname, 'AU3');
        XL_thous_ln_alpha   = xlswrite(XL_filename, export_thous_ln_alpha,     export_tabname, 'AV3');
    end
end
