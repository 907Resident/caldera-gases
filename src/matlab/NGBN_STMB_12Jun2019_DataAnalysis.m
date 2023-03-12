
%% NGBN_STMB_DataAnalysis

% NGBN_STMB_DataAnalysis
% Geoscientists: Ajayi, M. & Ayers, J.
% Location: Norris Geyser Basin - Steamboat (NGBN), Yellowstone National
% Park, WY
% GPS Coordinates: [44.723408, -110.703405]
% Data Type: Time Series Measurements 
% Metadata File: "...\Field Notes\Yellowstone\YNP_FieldNotes_06-15Jun2019.pdf"

%% Extra Functions Used in this Script
    % The following functions that were made by Ajayi, M.
    
        % lin_flux_AC()
        % KeelingCurve()
        % last10_iso
        % isocomp_array()
        % basic_gas_viz()
        % kp_viz()
        % Merge_PicData_Fxn()
        
%% Import Data Pre-Proccessed Data

% Create a directory to save files 
working_dir = "E:\moyoa\Documents\OneDrive - Vanderbilt\PhD_Dissertation\Data_Analysis\Picarro\Yellowstone\Jun2019\12Jun2019\";

% Date of measurements
ddmmmyyyy = "12Jun2019";

% Data should be imported via "Merge_PicData_Fxn"
 directory = "E:\moyoa\Documents\OneDrive - Vanderbilt\PhD_Dissertation\Data_Analysis\Picarro\Yellowstone\Jun2019\12Jun2019\*.dat";
[TT_PicData,PD_mtrx_nm,PD_mtrx_dt, dte] =                               ...
 Merge_PicData_Fxn(directory, 'America/Denver');  

%% Gather the NGBN_STMB Data from the Eosense Matrix 'PD_mtrx'

% Calculate the relative time for each measurement 
    % Pre-allocate vector for time normalized in milliseconds.
    % Milliseconds are used because of greater precision than seconds alone
    ms  = zeros([length(TT_PicData.Time), 1]);
    % Use for-loop to calculate the difference in times between each
    % measurement
    for i = 2:length(TT_PicData.Time)
            ms(i)           = ms(i-1) +                                 ...
            ( milliseconds(TT_PicData.Time(i) - TT_PicData.Time(i-1)) );
    end
    % Convert milliseconds back to seconds
    s                      = ms ./ 1E03;
    NGBN_STMB_rel_sec      = s;
    % Add the relative seconds vector to the timetable
    TT_PicData          = addvars(TT_PicData,NGBN_STMB_rel_sec,         ...
                         'Before','HP_CH4_dry');
    % Clear extraneous variables
    clearvars ms s

%% Spec Check
    % Ensure that measurements are not outside of the Picarro specs (e.g.
    % [CH4] <= 10 ppm & [CO2] <= 2500 ppm)

% Filter out Erroneous Measurements
TT_PicData.HP_CH4_dry(TT_PicData.HP_CH4_dry >= 1.00E01)         = NaN;
TT_PicData.HR_CH4_dry(TT_PicData.HR_CH4_dry >= 1.00E01)         = NaN;
TT_PicData.HP_d13_CH4(TT_PicData.HP_d13_CH4 >= 1.00E06)         = NaN;
TT_PicData.CO2_dry(TT_PicData.CO2_dry       >= 2.50E03)         = NaN;
TT_PicData.d13_CO2(TT_PicData.d13_CO2       >= 1.00E06)         = NaN;                             

%% Incorporate User Options

writeXL = input('Enter ''1'' to write main variables to an Excel file, enter zero if no Excel files are to be wrtitten: ');
% State how many times the eosAC cycled (opened and closed) during the
% long-term measurement
  % This will help to separate each closed chamber period during
  % that measurement session
nchams  = input('Enter the number of points measured at this site (reponses must be greater than zero): ');
    if nchams <= 0
       error('Please select a value greater than zero')
    elseif mod(nchams,1) ~= 0
        error('Please select an interger value')
    end
% Request a Site Tag to distinguish location from others
site_tag = input('Enter the four letter (all CAPS) code for the site you are analyzing: ', "s");
    % Ensure that the string is all upper case
    site_tag = upper(site_tag);
    
%%  Get Rid of the Duplicate Measurement Values
    % In some instances, users of the Picarro have observed duplicate
    % measurements. That is, measurements over 1-3 seconds that are exactly
    % the same.  Thus, it becomes necessary to remove these duplicate
    % measurements.

% The removal of the duplicate measurements can be accomplished by creating
% a vector 1's and 0's.  The 1's will be included in rows that are kept and
% the 0's will indicate rows that are to be removed.
    % For the NGBN_STMB, the data presented groups of four thus, only
    % every fourth row is retainted
RN = zeros([length(TT_PicData.HP_CH4_dry), 1]);
    % Assigns '1' to every fourth row in the vector RN
for i = 1:4:length(RN)
RN(i) = 1;
end
    % Create a new matrix that will contain the desired measurements and
    % the keep/delete indicator by horizontally concatenating RN to the
    % measurement vectors
RN_CH4 = [RN, datenum(TT_PicData.Time), NGBN_STMB_rel_sec               ...
         TT_PicData.HP_CH4_dry,         TT_PicData.HR_CH4_dry,          ...
         TT_PicData.HP_d13_CH4,         TT_PicData.HR_d13_CH4           ... 
         TT_PicData.CO2_dry,            TT_PicData.d13_CO2,             ...
         TT_PicData.Alarm_Status];
    % Retain every fourth entry and replace other entries with NaNs
for i = 1:length(RN)
    if RN_CH4(i,1) == 0
       RN_CH4(i,:) = NaN;
    end
end
    % Delete the rows that contain NaNs (as determined by the indicator
    % vector)
RN_CH4(any(isnan(RN_CH4),2),:)  = [];
    % Create a new vector that converts the edited time serial numbers to
    % datetime objects
Time                            = datestr(RN_CH4(:,2), -1);
Time                            = datetime(Time);

% Reassign the variables so that they represent the filtered information
    % Set up variables names
    VariableNames = {'Rel_Sec', 'HP_CH4_dry', 'HR_CH4_dry', 'HP_d13_CH4',  ...
                     'HR_d13_CH4', 'CO2_dry', 'd13_CO2', 'Alarm_Status'};
    % Recreate Timetable with pared down data, now called
    % *TT_NGBN_STMB*
    TT_NGBN_STMB = timetable(Time, RN_CH4(:,03),  RN_CH4(:,04),         ...
                           RN_CH4(:,05), RN_CH4(:,06), RN_CH4(:,07),    ...
                           RN_CH4(:,08), RN_CH4(:,09), RN_CH4(:,10),    ...
                          'VariableNames', VariableNames);
    % Record the summary of the data (includes descriptive statistics)
    smry_NGBN_STMB = summary(TT_NGBN_STMB);
    % Clear extraneous variables
    clearvars RN RN_CH4 RN_Time

%------Figure 01-------%
fig01 = figure;  
%------Figure 1A-------%
subplot(1,2,1)
    plot(TT_NGBN_STMB.Time, TT_NGBN_STMB.HP_CH4_dry, 'o')
        title('NGBN STMB 12Jun2019 [CH_4] vs. Timestamp (ALL)')
        ylabel('[CH_4] (ppm)')
        grid on
%         legend('[CH_4]', '± 1 std', 'Mean', 'Max', 'Min', ...
%            'Location', 'Northwest')

%------Figure 1B-------%
subplot(1,2,2)
    plot(TT_NGBN_STMB.Time, TT_NGBN_STMB.CO2_dry, 'mo')
        title('NGBN STMB 12Jun2019 [CO_2] vs. Timestamp (ALL)')
        ylabel('[CO_2] (ppm)')
        grid on
%     legend('[CO_2]', '± 1 std', 'Mean', 'Max', 'Min', ...
%            'Location', 'Northwest')
% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_timeline_side-by-side.fig", ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig01, fig_file)

%------Figure 02-------%
fig02 = figure;
left_ax_color          = [0 0.477 0.741];
right_ax_color         = [1 0 1];

set(fig02,'defaultAxesColorOrder',[left_ax_color; right_ax_color]);

yyaxis left
% --- CH4 --- %
   h1 = plot(TT_NGBN_STMB.Time,TT_NGBN_STMB.HP_CH4_dry, 'o');
                    title('NGBN STMB 12Jun2019 [CH_4] and [CO_2] vs. Time')
                    grid on
                    % Modify the line markers
                    h1.Marker               = 'o';
                    h1.LineStyle            = 'none';
                    ylabel('[CH_4] (ppm)')
yyaxis right                    
% --- CO2 --- %
   h2 = plot(TT_NGBN_STMB.Time,TT_NGBN_STMB.CO2_dry, 'o');
                    h2.Marker               = 'o';
                    h2.LineStyle            = 'none';
                    h2.Color                = [1 0 1];
                    ylabel('[CO_2] (ppm)')

                    % Add a legend
                    legend({'[CH_4]', '[CO_2]'}, 'EdgeColor', 'none',   ...
                            'Color','none', 'Location', 'Best')
% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_timeline_overlay.fig",  ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig02, fig_file)  
                        
% Clear extraneous variables
clearvars h1 h2 left_ax_color right_ax_color
clearvars PicData PD_mtrx_dt TT_PicData 

%% Extract and Isolate Pre/ChamON/Post Data
    % Here, we separate the data collected into packages for before,
    % during, and after chamber enclosure.  This requires visual inspection
    % of the CO2 data (see first two graphs) and comparison to the field
    % notes associated with this measurement session

% Create array of time contraints 
DATE                            = dte(1:nchams);

%----PreBack----%
    % Pre-allocate the datetime array
PreBack_Times                   = NaT([nchams, 2]);

PreBack_Times(:,1)              = datetime(                             ...
['10:13:19'; '10:35:00'; '10:54:58'; '11:14:52'; '11:34:56'; '11:55:02';...
 '12:14:58'; '12:34:45'; '12:54:52'; '13:11:46'; '13:32:10'; '13:51:43';...
 '14:11:54'; '14:31:51'; '14:51:42'; '15:11:42'; '15:32:14'],           ...
 'InputFormat', 'HH:mm:ss');
    % Extract the time from the combined date-time value
PreBack_Times(:,1)              = DATE + timeofday(PreBack_Times(:,1));

PreBack_Times(:,2)              = datetime(                             ...
['10:18:06'; '10:39:22'; '10:59:16'; '11:18:38'; '11:38:56'; '11:58:49';...
 '12:18:56'; '12:39:20'; '12:56:09'; '13:15:47'; '13:35:57'; '13:55:09';...
 '14:15:23'; '14:35:28'; '14:55:19'; '15:15:38'; '15:35:53'],           ...
 'InputFormat', 'HH:mm:ss');
    % Extract the time from the combined date-time value
PreBack_Times(:,2)              = DATE + timeofday(PreBack_Times(:,2));

%----ChamON----%
ChamON_Times                    = NaT([nchams, 2]);

ChamON_Times(:,1)               = datetime(                             ...
['10:18:10'; '10:39:26'; '10:59:19'; '11:18:41'; '11:39:00'; '11:58:53';...
 '12:19:00'; '12:39:24'; '12:56:12'; '13:15:50'; '13:36:00'; '13:55:22';...
 '14:15:26'; '14:35:31'; '14:55:22'; '15:15:42'; '15:35:57'],           ...
 'InputFormat', 'HH:mm:ss');

ChamON_Times(:,1)               = DATE + timeofday(ChamON_Times(:,1));

% Fill in with the **"end"** time of each chamber closure period (i.e. the
% time when the chamber opens
ChamON_Times(:,2)               = datetime(                             ...
['10:34:03'; '10:54:03'; '11:14:00'; '11:34:00'; '11:53:59'; '12:13:53';...
 '12:33:54'; '12:53:57'; '13:10:47'; '13:30:56'; '13:50:47'; '14:10:51';...
 '14:30:45'; '14:50:43'; '15:10:42'; '15:30:54'; '15:50:46'],           ...
 'InputFormat', 'HH:mm:ss');

ChamON_Times(:,2)               = DATE + timeofday(ChamON_Times(:,2));

%----PostBack----%
    % No post background data needed here; the pre background is used as
    % pre and post for these continuous measurements
    
    % Clear extraneous variables
    clearvars DATE
    
%% Extract the Individual Measurement
    % Use the start and end time from the previous section, use the
    % TIMERANGE() function to isolate the data.  Then transfer the results
    % from a timetable to a multidimensional array

    NGBN_STMB_PreBack    = NaN([100, width(TT_NGBN_STMB)+1, nchams]);
    NGBN_STMB_ChamON     = NaN([420, width(TT_NGBN_STMB)+1, nchams]);
    
for i = 1:nchams
    % PreBack
    S_Pre           = timerange(PreBack_Times(i,1), PreBack_Times(i,2), ...
                               'closed');
    T               = TT_NGBN_STMB(S_Pre,:);
    A               = NaN([height(T), width(T)+1]);
    dnum            = datenum(T.Time);
    A(:,1)          = dnum;
    A(:,2:end)      = table2array(T);
   [m, n]           = size(A);
    
    NGBN_STMB_PreBack(1:m,1:n,i) = A;
    
    % For posterity clear extraneous variables here, not required but in
    % case future changes to the code are made, this line will eliminate
    % mix ups
    clearvars T A dnum m n
    
    % ChamON
    S_ChamON    = timerange(ChamON_Times(i,1), ChamON_Times(i,2),       ...
                            'closed');
    T           = TT_NGBN_STMB(S_ChamON,:);
    A           = NaN([height(T), width(T)+1]);

    dnum        = datenum(T.Time);
    A(:,1)      = dnum;
    A(:,2:end)  = table2array(T);
    % Change the seconds to relative seconds of each individual measurement
    % (0 to ~1000/1200)
        % Calculate the differences in time (sec) between each measurement
        dmy_diff= diff(A(:,2));
        % Usimg CUMSUM(), calculate the accumlated time for each
        % measurement.  In other words, this will now show the relative
        % within a measurement period (e.g. the measurement occured 45 sec
        % into the measuring period)
        A(2:end,2)  = cumsum(dmy_diff);
        % Set the first time component to zero because no time has elapsed
        % at the first measurement within the time peropd
        A(1,2)  = 0;
    % Acquire the dimnensions of the array A so that it can be properly
    % placed in the ___ data
   [m, n]       = size(A);
    
    NGBN_STMB_ChamON(1:m,1:n,i) = A;  
    
    % Clear extra variables
    clearvars T A dnum m n dmy_diff
    
end

%% Extract Datetimes

    % PreBack
        % Pre-allocate
    NGBN_STMB_PreBack_DateTimes     = NaT([100 1 nchams]);
for idx = 1:nchams
    f                               = find(isnan(NGBN_STMB_PreBack(:,1,idx)));
    NGBN_STMB_PreBack_DateTimes(1:min(f)-1,1,idx)    = datetime(        ...
         datestr(NGBN_STMB_PreBack(1:min(f)-1,1,idx)));
end

    % ChamON
        % Pre-allocate
    NGBN_STMB_ChamON_DateTimes      = NaT([350 1 nchams]);
for idx = 1:nchams
    f                               = find(isnan(NGBN_STMB_ChamON(:,1,idx)));
    NGBN_STMB_ChamON_DateTimes(1:min(f)-1,1,idx)    = datetime(         ...
         datestr(NGBN_STMB_ChamON(1:min(f)-1,1,idx)));
end

% Delete Extraneous Variables
clearvars S_Pre S_ChamON

%% Spec Check
    % Ensure that measurements are not outside of the Picarro specs (e.g.
    % [CH4] <= 10 ppm & [CO2] <= 2500 ppm)

% CH4 Spec Check
for sp_idx = 1:nchams
    CH4_Sp_Chk = find(NGBN_STMB_ChamON(:,3,sp_idx) > 10);
    if  any(CH4_Sp_Chk)
        NGBN_STMB_ChamON(CH4_Sp_Chk,3,sp_idx) = NaN;
        NGBN_STMB_ChamON(CH4_Sp_Chk,5,sp_idx) = NaN;
    else
        Spec_Chk_str = sprintf('There were no [CH4] measurements greater than 10 ppm for measurement number %d.', ...
                                                                   sp_idx);
        disp(Spec_Chk_str)                                                       
    end
end

% CO2 Spec Check
for sp_idx = 1:nchams
    CO2_Sp_Chk = find(NGBN_STMB_ChamON(:,7,sp_idx) > 2500);
    if  any(CO2_Sp_Chk)
        NGBN_STMB_ChamON(CO2_Sp_Chk,7:8,sp_idx) = NaN;
    else
        Spec_Chk_str = sprintf('There were no [CO2] measurements greater than 2500 ppm for measurement number %d.', ...
                                                                   sp_idx);
        disp(Spec_Chk_str)                                                       
    end
end

clearvars CH4_Sp_Chk CO2_Sp_Chk sp_idx

%% Prelinary Visualizations

[fig03, fig04, fig05] = basic_gas_viz(NGBN_STMB_ChamON, nchams, site_tag,    ...
                               ddmmmyyyy, working_dir);

% Save figure 03 as .fig to working directory
fi       = sprintf("%s_%s_CH4_and_CO2_timeline_array.fig",              ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig03, fig_file) 
% Save Figures 04-05 as .jpg manually or as needed
                           
% Clear extraneous variables
clearvars CH4_pl CO2_pl title_str pnt trans

%% Convert gas measurements from ppm tg cgs units

% Converting constants
Vm                              = 22.71108;  % Standard Molar Volume (mol L^-1)
CH4_molec_wt                    = 16.043;    % g/mol
CO2_molec_wt                    = 44.009;    % g/mol
% Imported concentrations (resampled)
c_CH4_ppm                       = NGBN_STMB_ChamON(:,3,:);   
c_CO2_ppm                       = NGBN_STMB_ChamON(:,7,:);
% Chamber Height (m)
    % 5in = 12.70 cm = 0.127 m
H                               = 0.127; 
% Number of observations (N)
N                               = length(c_CH4_ppm);

%%% To get fluxes, think of concentration expressed in micrograms per cubic
%%% meter ug/m3
    %%% Convert from ratio of gas to air (ppm) to mass per volume (mg /
    %%% m^3)
        %%% (Molec_wt / Vm) * (1000 L / m^3) * ppm = ug / m^3
        c_CH4_cgs               = ((CH4_molec_wt / Vm) .* 1E03)         ...
                                    .* c_CH4_ppm;
        %%%% (ug / m^3) / 1000 = mg / m^3 
        c_CH4_cgs               = c_CH4_cgs ./ 1E03;
    %%% Repeat for carbon dioxide     
        %%% (Molec_wt / Vm) * (1000 L / m^3) * ppm = ug / m^3
        c_CO2_cgs               = ((CO2_molec_wt / Vm)  .* 1E03)        ...
                                    .* c_CO2_ppm;
        %%%% (ug / m^3) / 1000 =  mg / m^3 
        c_CO2_cgs               = c_CO2_cgs ./ 1E03;        

% Add the converted cgs concentrations to the main numeric array
 for idx = 1:nchams
        NGBN_STMB_ChamON(:,10,idx)    = c_CH4_cgs(:,1,idx);
        NGBN_STMB_ChamON(:,11,idx)    = c_CO2_cgs(:,1,idx);
 end

% Clear extraneous variables
clearvars N Vm CH4_molec_wt CO2_molec_wt c_CH4_cgs_UnEx c_CO2_cgs_UnEx c_CH4_ppm c_CO2_ppm 

%% Flux Calculation

[lin_flux, lin_mdl, lin_slope, flx_srt_time, flx_end_time, fig06, fig07] = ...
 lin_flux_AC(NGBN_STMB_ChamON, "y");
% Save figure as .fig to working directory
fi       = sprintf("MATLAB_Figs\\%s_%s_CH4_lin_mdl_array.fig",                       ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig06, fig_file)

% Save figure as .fig to working directory
fi       = sprintf("MATLAB_Figs\\%s_%s_CO2_lin_mdl_array.fig",                       ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig07, fig_file)

% Save figures as .jpg to working directory
fi       = sprintf("%s_%s_CH4_lin_mdl_array.jpg",                       ...
                    site_tag, ddmmmyyyy);
jpg_file = working_dir+fi;
saveas(fig06, jpg_file)

% Save figure as .jpg to working directory
fi       = sprintf("%s_%s_CO2_lin_mdl_array.jpg",                       ...
                    site_tag, ddmmmyyyy);
jpg_file = working_dir+fi;
saveas(fig07, jpg_file)
    
%----Exponential Model Analysis----%
    % A two-minute equilibration period (approx 33-36 measurements) is
    % allowed before an assessment of the exponential fit.  After this
    % two-minute period, the remaining portion of the measurement period is
    % evaluted using an exponential fit:
        % C_t = psi + (C_0 - psi) * exp(-kappa*t)
    % Pre-allocate vector for fluxes (per measurement cycle)
    exp_mdl     = cell([1 nchams 2]);
    % Define the options for the exponential fit (i.e. fitting method,
    % starting point for the iteration)
    exp_ft_ops  = fitoptions('StartPoint', [0 0],                       ...
                             'Method', 'NonlinearLeastSquares');
    % Create an exponential model to evaluate the flux data using fittype()
    exp_fxn     = fittype(@(psi, kappa, c_0, t) psi + (c_0 - psi)       ...
                          * exp(-kappa*t),                              ...
                          'dependent', {'C_t'}, 'independent', {'t'},   ...
                          'coefficients', {'psi', 'kappa'}, 'problem',  ...
                          'c_0', 'options', exp_ft_ops);
        % (This is equivalent to 1.9 ppm CH4 and could be adjusted to
        % indicate the concentration at time = 0)
    exp_slope   = zeros([1  3       nchams      2]);
    exp_flux    = zeros([1  nchams              2]);
    
    x_expData   = cell([1   nchams  2            ]);
    y_expData   = cell([1   nchams  2            ]);

for i = 1:nchams
%---CH4---%
% Apply an exponential model to evaluate the flux data using fit()
    % The row indicies are hard coded to ensure that no NaN values are
    % analyzed; additionally, this segment of data accounts for the
    % majority of the data after the equilibriation period
[x_exp_prepData, y_exp_prepData]    =                                   ...
                 prepareCurveData(NGBN_STMB_ChamON(:,02,i),             ...
                                  NGBN_STMB_ChamON(:,10,i));   
x_expData{:, i, 1}                  = x_exp_prepData;
y_expData{:, i, 1}                  = y_exp_prepData;
exp_mdl  {1, i, 1}                  = fit(x_expData{i},y_expData{i},    ...
                                          exp_fxn,'problem',            ...
                                          NGBN_STMB_ChamON(01,10,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
exp_slope(1, 1, i, 1)               = exp_mdl{i}.psi;
exp_slope(1, 2, i, 1)               = exp_mdl{i}.kappa;
exp_slope(1, 3, i, 1)               = exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
exp_flux (1, i, 1)  =((exp_mdl{1, i, 1}.psi - exp_mdl{1, i, 1}.c_0) .*  ...
                       exp_mdl{1, i, 1}.kappa) .* H .* 3600;
                   
%---CO2---%
[x_exp_prepData, y_exp_prepData]    =                                   ...
                 prepareCurveData(NGBN_STMB_ChamON(:,02,i),             ...
                                  NGBN_STMB_ChamON(:,11,i));   
x_expData{:, i, 2}                  = x_exp_prepData;
y_expData{:, i, 2}                  = y_exp_prepData;
exp_mdl  {1, i, 2}                  = fit(x_expData{i},y_expData{i},    ...
                                          exp_fxn,'problem',            ...
                                          NGBN_STMB_ChamON(01,11,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
exp_slope(1, 1, i, 2)               = exp_mdl{i}.psi;
exp_slope(1, 2, i, 2)               = exp_mdl{i}.kappa;
exp_slope(1, 3, i, 2)               = exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
exp_flux (1, i, 2)  =((exp_mdl{1, i, 2}.psi - exp_mdl{1, i, 2}.c_0) .*  ...
                       exp_mdl{1, i, 2}.kappa) .* H .* 3600;
end

% Reshape Exponential Flux to Match Shape of Linear Flux matrix [m x 1]
exp_flux = reshape(exp_flux,[length(exp_flux) 1, 2]);

% Visualize the fluxes
fig08 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    %  CH4
    nexttile
        % Error Sactter
        eb_flux = errorbar(1:nchams, lin_flux(:,1,1),                   ...
                           lin_flux(:,1,1) - lin_flux(:,2,1),           ...
                           lin_flux(:,3,1) - lin_flux(:,1,1),           ...
                           'o', 'CapSize', 0, 'LineWidth', 1.25,        ...
                           'Color', 'k'); 
        set(eb_flux, 'MarkerFaceColor', "#4da6ff", 'MarkerEdgeColor', 'k'),...
            grid on
%         if any(abs(exp_flux(:,:,1) - lin_flux(:,:,1)) <= 1E02) == 1
%             hold on
%             bar(1:nchams,exp_flux(:,:,1), 0.5, 'FaceColor', 'c')
%         else
%             disp('The absolute difference between the exponential and linear models was greater than 100 mg m^-2 hr^-1 for all measuerments of CH4. Reexamine the models.')
%         end
        tstr = sprintf("%s %s | CH_4 Fluxes", site_tag, ddmmmyyyy);
        title(tstr)
        xlabel('Location')
        ylabel('F_{CH_4} (mg m^{-2} hr^{-1})')
%       legend({'Linear Model', 'Exponential Model'}, 'Location', 'SW')
        xticks(1:nchams)
        xticklabels({'1.1', '1.2', '1.3', '1,4', '1.5', '1.6', '1.7', '1.8',...
                     '1.9', '1.10', '1.11', '1.12', '1.13', '1.14', '1.15', ...
                     '1.16', '1.17'})
        hold off
        % Density Plot
        nexttile
       [f, xi] = ksdensity(lin_flux(:,1,1));
       plot(xi,f, 'MarkerFaceColor', "#4da6ff",                         ...
                            'MarkerEdgeColor', 'k');
       xlabel('F_{CH_4} (mg m^{-2} hr^{-1})'), ylabel("Density")
       title(tstr)
       grid on
% CO2
    nexttile
        % Error Sactter
        eb_flux = errorbar(1:nchams, lin_flux(:,1,2),                   ...
                           lin_flux(:,1,2) - lin_flux(:,2,2),           ...
                           lin_flux(:,3,2) - lin_flux(:,1,2),           ...
                           'o', 'CapSize', 0, 'LineWidth', 1.25,        ...
                           'Color', 'k'); 
        set(eb_flux, 'MarkerFaceColor', "#ff66ff", 'MarkerEdgeColor', 'k'),...
            grid on
%         if any(abs(exp_flux(:,:,2) - lin_flux(:,:,2)) <= 1E02) == 1
%             hold on
%             bar(1:nchams,exp_flux(:,:,2), 0.5,                      ...
%                 'FaceColor', [102/255 0/255 102/255])
%         else
%             disp('The absolute difference between the exponential and linear models was greater than 100 mg m^-2 hr^-1 for all measurements of CO2. Reexamine the models.')
%         end
        tstr = sprintf("%s %s | CO_2 Fluxes", site_tag, ddmmmyyyy);
        title(tstr)
        xlabel('Location')
        ylabel('F_{CO_2} (mg m^{-2} hr^{-1})')
        xticks(1:nchams)
        xticklabels({'1.1', '1.2', '1.3', '1,4', '1.5', '1.6', '1.7', '1.8',...
                     '1.9', '1.10', '1.11', '1.12', '1.13', '1.14', '1.15', ...
                     '1.16', '1.17'})
        % Density Plot
        nexttile
       [f, xi] = ksdensity(lin_flux(:,1,2));
       plot(xi,f, 'm', 'MarkerFaceColor', "#ff66ff",                    ...
                       'MarkerEdgeColor', 'k');
       xlabel('F_{CO_2} (mg m^{-2} hr^{-1})'), ylabel("Density")
       title(tstr)
       grid on        
% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_fluxes_bar.fig",     ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;                
savefig(fig08, fig_file)

% Save figure as .jpg to working directory
fi       = sprintf("%s_%s_CH4_and_CO2_fluxes_scatter_and_density.jpg",  ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;                
saveas(fig08, fig_file)       

%% Keeling Plots
    % Keeling plots (originally generated by Keeling (1961(?)) to assess
    % the source of the gas

% Pre-allocate vectors for KeelingCurve Results
KP_rslts                        = cell( [1 nchams   2        ]);
KP_coeffs                       = zeros([1 2        nchams  2]);
KP_coeffs_ci                    = zeros([2 2        nchams  2]);
KP_Fits                         = cell([nchams 4]);

for idx = 1:nchams
%--CH4--%
    % Prepare data for the regression (i.e. remove NaN values, shape the
    % vectors as vertical vectors, etc.) 
[xData, yData]  = prepareCurveData((1 ./ NGBN_STMB_ChamON(:,3,idx)),    ...
                                         NGBN_STMB_ChamON(:,5,idx));
    % Set up fittype and options.
KP_ft = fittype('poly1');

    % Fit model to data.
[KP_fitresult,KP_gof]           = fit(xData, yData, KP_ft);
 KP_Fits{idx,1}                 = KP_fitresult;
 KP_Fits{idx,2}                 = KP_gof.rsquare;

    % Gather the coefficients for the fit equation
fitvals = coeffvalues(KP_fitresult);
    % Gather the confidence intervals for c
    KP_rslts{1,idx,1}           = KP_fitresult;
    KP_coeffs(1,1,idx,1)        = fitvals(1);
    KP_coeffs(1,2,idx,1)        = fitvals(2);
    KP_coeffs_ci(:,:,idx,1)     = confint(KP_fitresult, 0.95);
    
%--CO2--%
[xData, yData]  = prepareCurveData((1 ./ NGBN_STMB_ChamON(:,7,idx)),    ...
                                         NGBN_STMB_ChamON(:,8,idx));
    % Set up fittype and options.
KP_ft                           = fittype('poly1');

    % Fit model to data.
[KP_fitresult,KP_gof]           = fit(xData, yData, KP_ft);
 KP_Fits{idx,3}                 = KP_fitresult;
 KP_Fits{idx,4}                 = KP_gof.rsquare;

    % Gather the coefficients for the fit equation
fitvals                         = coeffvalues(KP_fitresult);
    % Gather the confidence intervals for c
    KP_rslts{1,idx,1}           = KP_fitresult;
    KP_coeffs(1,1,idx,2)        = fitvals(1);
    KP_coeffs(1,2,idx,2)        = fitvals(2);
    KP_coeffs_ci(:,:,idx,2)     = confint(KP_fitresult, 0.95);
end

% Visualize the Keeling plots
[fig09, fig10] = kp_viz(NGBN_STMB_ChamON, nchams,                       ...
                        KP_Fits, KP_coeffs, KP_coeffs_ci, site_tag,     ...
                        ddmmmyyyy, working_dir);

%% Isotope Calculations using "Last 10" Method

% Calculate the isotopic composition using the "Last 10" method as
% described in the help content for last10_iso()
[last10_d13CH4,  last10_d13CO2, fig11] = last10_iso(NGBN_STMB_ChamON, "y");

% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_d13C_last10_EDA.fig",...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi; 
savefig(fig_file)                    

% Save figure as .jpg to working directory
fi       = sprintf("%s_%s_CH4_and_CO2_d13C_last10_EDA.jpg",             ...
                    site_tag, ddmmmyyyy);
jpg_file = working_dir+fi;
saveas(fig11, jpg_file)

%% Calculate alpha
    % a = [(?13C-CO2 - ?13C-CH4)/1000 ] + 1
    Alpha           = zeros([1,nchams]);
    thous_ln_alpha  = zeros([1,nchams]);
    for i = 1:nchams
        Alpha(i)            = ((KP_coeffs(1,2,i,2) - KP_coeffs(1,2,i,1))...
                                ./ 1E03) + 1;
        thous_ln_alpha(i)   = 1E03 .* log(Alpha(i));
    end

    % Calculate alpha using the method in Forde et al. 2019
        % a = [1E03 + ?13C-CO2] / [1E03 + ?13C-CH4]
    alpha_Forde19   = zeros([1 nchams]);
    for i = 1:nchams
        alpha_Forde19(i)    = (1E03 + KP_coeffs(1,2,i,2)) ./            ...
                              (1E03 + KP_coeffs(1,2,i,1) );
    end

%% Carbon Isotope Comparisons
    % Using a chart (originally generated by Whiticar (1999) we can further
    % elucidate the source of the gas via the carbon isotopes of the gas

% Isolate the isotopic data
    iCH4 = NGBN_STMB_ChamON(:,5,:); iCO2 = NGBN_STMB_ChamON(:,8,:);
    % Use the isocomp_array() to plot the iCO2 vs. iCH4 data
    [~, ~, iso_array]   = isocomp_array(iCH4, iCO2, nchams, "Zoom", "yes");
    
    % Save figure as .fig to working directory
    fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_d13C_isocomp.fig",  ...
                        site_tag, ddmmmyyyy);
    fig_file = working_dir+fi; 
    savefig(iso_array, fig_file)    
    
    % Save figure as .jpg to working directory
    fi       = sprintf("%s_%s_CH4_and_CO2_d13C_isocomp.jpg",            ...
                        site_tag, ddmmmyyyy);
    jpg_file = working_dir+fi;
    saveas(iso_array, jpg_file)
    
%% CO2/CH4 Ratios

% Create an array of CH4/CO2 vs Time for each point measured
NGBN_STMB_ChamON_Ratio           = NGBN_STMB_ChamON(:,7,:) ./ NGBN_STMB_ChamON(:,3,:);
% Record the descriptive stats for the CH4/CO2 ratios
    % ... _stats = [mean; median; range; standard deviation];
    % *** COMPLETE the block BELOW (21 Aug 2018) ***
NGBN_STMB_ChamON_Ratio_stats     = [mean(NGBN_STMB_ChamON_Ratio);       ...
                                    median(NGBN_STMB_ChamON_Ratio);     ...
                                    range(NGBN_STMB_ChamON_Ratio);      ...
                                    std(NGBN_STMB_ChamON_Ratio)];

    % CO2/CH4
fig12 = figure;
for idx = 1:nchams

    nexttile
        s_CH4 = scatter(NGBN_STMB_ChamON(:,2,idx), NGBN_STMB_ChamON_Ratio(:,idx), ...
                        36, [127/255 0/255 255/255]);
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
            if idx     <= nchams
                trans   = 1;
                pnt     = idx;
            else
                trans   = 99;
                pnt     = 99;
            end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str   = sprintf(                                  ...
                             'NGBN STMB %d.%d %s [CO_2]/[CH_4] vs. Rel Time',  ...
                              trans,pnt,ddmmmyyyy);
                title(title_str, 'FontSize', 8)
                % The maximum value on the x-axis will be 1200 (= 60 sec *
                % 20 min)
                xlim([0 1200])
                % Provide tick marks in 300 sec (5 min) intervals
                xticks(0:300:1200)
                % Add labels to axis
                xlabel('Relative Time (sec)', 'FontSize', 10)
                ylabel('[CO_2]/[CH_4]', 'FontSize', 10)
                % Apply gridding on both axis to see the data better
                grid on
end

    % Save figure as .fig to working directory
    fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_ratio_array.fig",   ...
                        site_tag, ddmmmyyyy);
    fig_file = working_dir+fi;
    savefig(fig12, fig_file)
%---- Scatter plots of carbon isotopes and CO2-CH4 Ratio ----%

fig13 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for i = 1:nchams
    
     % Get the for each chamber placement
     X = NGBN_STMB_ChamON(:,5,i); Y = NGBN_STMB_ChamON(:,8,i); 
     % Color the data according to the CO2-CH4 ratio
     color = NGBN_STMB_ChamON_Ratio(:,i);
     
     % Create a grid of plots for each chamber measurement
     nexttile
     % Use the SCATTER() function
     s13 = scatter(X,Y,[],color, 'filled');
     set(s13, "MarkerEdgeColor", "k")
     grid on
     % This if-elseif-else block is used to appropriately designate the
     % numerical portion of each title in the array of plots
            if idx     <= nchams
                trans   = 1;
                pnt     = idx;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str   = sprintf(                                         ...
            '%s %d.%d %s d^{13}C-CO_2 vs. d^{13}C-CH_4',                ...
             site_tag,trans,pnt,ddmmmyyyy);
             title(title_str, 'FontSize', 8)
     % Use the "Cool" colormap
     colormap cool
     c = colorbar;
        c.Label.String = '[CO_2]/[CH_4]';
       %c.Limits       = [min(NGBN_STMB_ChamON_Ratio) max(NGBN_STMB_ChamON_Ratio)];
     
end

% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4-CO2_d13C_ratio_array.fig",   ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig13, fig_file)

% Save figure as .jpg to working directory
fi       = sprintf("%s_%s_CH4-CO2_d13C_ratio_array.jpg",                ...
                    site_tag, ddmmmyyyy);
jpg_file = working_dir+fi;
saveas(fig13, jpg_file)

%Clear extraneous variables
clearvars title_str trans pnt X Y color
    
%% 3-Dimensional Analysis

% Create a plot matrix of 3-D plots using d13C-CH4 (x-axis), d13C-CO2
% (y-axis), [CO2]-[CH4] ratio (z-axis)

    % Organize the data into the vertical vectors for the x-, y-, and
    % z-axis.  This will be a matrix with with a length of  , three columns
    % representing the metrics corresponding to the three axis, and a
    % third dimension in which each "page" represents an individual
    % chamber measurement.  The layout of the "xyz" matrix is:
    % [num of obs, x/y/z axis, cham meas (1:ncham)]
    xyz = NaN([length(NGBN_STMB_ChamON), 3, nchams]);
    for i = 1:nchams
        xyz(:,1,i) = NGBN_STMB_ChamON(:,5,i);
        xyz(:,2,i) = NGBN_STMB_ChamON(:,8,i);
        xyz(:,3,i) = NGBN_STMB_ChamON_Ratio(:,1,i);
    end

% Pre-allocate the surface fits for the 3-D plots
    sf3 = cell([nchams, 3]);
    % 3-D Plots
fig14 = figure;
for idx = 1:nchams
    % Extract Data for each chamber measurement
    data_3D         = xyz(:,:,idx);
    
    % Eliminate the NaNs 
    data_3D_no_NaNs = data_3D(~any(isnan(data_3D),2),:);
    
    % Create a simple polynomial surface in the X & Y directions for the
    % plot3 charts
    [sf3{idx,1}, sf3{idx,2}, sf3{idx,3}]                                ...
                    = fit([data_3D_no_NaNs(:,1), data_3D_no_NaNs(:,2)], ...
                           data_3D_no_NaNs(:,3), 'poly22');
    % Create subplot matrix
    nexttile
        p3_CH4              =  plot(sf3{idx},                           ...
                              [data_3D_no_NaNs(:,1), data_3D_no_NaNs(:,2)],...
                               data_3D_no_NaNs(:,3));
       % Apply gridding on both axis to see the data better
                       grid on
                       
                       set(gcf, 'color', [0.80 0.80 0.80]);
                       colormap  cool
                       alpha= 0.65;
                       view(-62.14,32.56)
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if idx     <= nchams
                    trans   = 1;
                    pnt     = idx;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str   = sprintf('NGBN STMB %d.%d %s Conc and Isotope Data',...
                                       trans, pnt, ddmmmyyyy);
                title(title_str, 'FontSize', 7)
                % Add labels to axis
                xlabel('\delta^{13}C-CH_{4} (‰)', 'FontSize', 8)
                ylabel('\delta^{13}C-CO_{2} (‰)', 'FontSize', 8)
                zlabel('[CO_2]/[CH_4]', 'FontSize', 8)
                
end

    % Save figure as .fig to working directory
    fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_ratio_array.fig",   ...
                        site_tag, ddmmmyyyy);
    fig_file = working_dir+fi; 
    savefig(fig14, fig_file)

% Plot residuals of the surface fit
fig15 = figure;
for i = 1:nchams
    
     nexttile
        hist_resid          = histogram(sf3{i,3}.residuals,             ...
                              'NumBins', 8,                             ...
                              'normalization', 'pdf');
    
    % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if i       <= nchams
                    trans   = 1;
                    pnt     = i;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str   = sprintf(                                  ...
                             'NGBN STMB %d.%d %s Residual of Surface Fit', ...
                              trans,pnt,ddmmmyyyy);
                title(title_str, 'FontSize', 7)
                % Add labels to axis
                xlabel('Residual', 'FontSize', 8)
                ylabel('Probability Density', 'FontSize', 7)
                grid on
end

% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4-CO2_ratio_hist_resids_array.fig",...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig15, fig_file)

%Clear extraneous variables
clearvars title_str trans pnt data_3D data_3D_no_NaNs

%% Perimeter EGM Measurement


%% Pre- and Post-Background EDA
    % This section provides some quick EDA for the pre- and post-background
    % (chamber open between measurements)
    
% Run basic stats on the pre and post data    
background_CH4 = cell([nchams 2]);      
background_CO2 = cell([nchams 2]);   

% Convert all values less than or equal to zero to NaN
NGBN_STMB_PreBack(NGBN_STMB_PreBack == 0)     = NaN;

for i = 1:nchams 
    % Descripitive Statistics
        % CH4
    background_CH4{i,1}          = desc_stats(NGBN_STMB_PreBack(:,3,i));
    background_CH4{i,2}          = desc_stats(NGBN_STMB_PreBack(:,5,i));
         % CO2
    background_CO2{i,1}          = desc_stats(NGBN_STMB_PreBack(:,7,i));
    background_CO2{i,2}          = desc_stats(NGBN_STMB_PreBack(:,7,i));   
    
    % Concentrations
    fig = figure;
        % [CH4]
    nexttile
        % Pre-Background
    [f, xi, bw_pre] = ksdensity(NGBN_STMB_PreBack(:,3,i));
    plot(xi, f)
    xlabel("[CH_{4}] (ppm)", "Interpreter", "tex")
    ylabel("Density")    
    % Add annotation if not running script to be exported MS Excel
    if writeXL ~= 1
            % Texboxt annotation
        tb_str = sprintf("Bandwidth Estimator for Pre-Back = %2.3f", bw_pre);
        dim    =[0.7848958,0.8506750,0.144271,0.045691];
        tb     = annotation('textbox', dim, 'String', tb_str, 'FitBoxToText','on');
        tb.EdgeColor = "none";
    end
    grid on

        % [CO2]
    nexttile
        % Pre-Background
    [f, xi, bw_pre] = ksdensity(NGBN_STMB_PreBack(:,7,i));
    plot(xi, f, 'm')
    xlabel("[CO_{2}] (ppm)", "Interpreter", "tex")
    ylabel("Density")    
    % Add annotation if not running script to be exported MS Excel
    if writeXL ~= 1
            % Texboxt annotation
        tb_str = sprintf("Bandwidth Estimator for Pre-Back = %2.3f", bw_pre);
        dim    =[0.7848958,0.8506750,0.144271,0.045691];
        tb     = annotation('textbox', dim, 'String', tb_str, 'FitBoxToText','on');
        tb.EdgeColor = "none";
    end
    grid on

        % d13CH4
    nexttile
        % Pre-Background
    [f, xi, bw_pre] = ksdensity(NGBN_STMB_PreBack(:,5,i));
    plot(xi, f)
    xlabel("δ^{13}C-CH_{4} (‰)", "Interpreter", "tex")
    ylabel("Density")    
    % Add annotation if not running script to be exported MS Excel
    if writeXL ~= 1
            % Texboxt annotation
        tb_str = sprintf("Bandwidth Estimator for Pre-Back = %2.3f", bw_pre);
        dim    =[0.7848958,0.8506750,0.144271,0.045691];
        tb     = annotation('textbox', dim, 'String', tb_str, 'FitBoxToText','on');
        tb.EdgeColor = "none";
    end
    grid on 

        % d13CO2
    nexttile
        % Pre-Background
    [f, xi, bw_pre] = ksdensity(NGBN_STMB_PreBack(:,8,i));
    plot(xi, f, 'm')
    xlabel("δ^{13}C-CO_{2} (‰)", "Interpreter", "tex")
    ylabel("Density")       
    % Add annotation if not running script to be exported MS Excel
    if writeXL ~= 1
            % Texboxt annotation
        tb_str = sprintf("Bandwidth Estimator for Pre-Back = %2.3f", bw_pre);
        dim    =[0.7848958,0.8506750,0.144271,0.045691];
        tb     = annotation('textbox', dim, 'String', tb_str, 'FitBoxToText','on');
        tb.EdgeColor = "none";
    end
    grid on
    
    % Save figure as .fig to working directory
    file_str = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_background_EDA_ncham-%d.fig", ...
                        site_tag, ddmmmyyyy, i);
    fig_file = working_dir+file_str;
    savefig(fig, fig_file)
    
    % Save figure as .jpg to working directory
    file_str = sprintf("%s_%s_CH4_and_CO2_background_EDA_ncham-%d.jpg", ...
                        site_tag, ddmmmyyyy, i);
    fig_file = working_dir+file_str;
    saveas(fig, fig_file)
    
end

%% Export Data
    
% Write the pre-background, chamber, and post-background data into an Excel
% file to be used by non-MATLAB users who are a part of the project

if writeXL == 1
    for i = 1:nchams
        % Create a string to label the  tab of the excel sheet
            if i               <= nchams
                trans           = 1;
                pnt             = i;
            else
                trans           = 99;
                pnt             = 99;
            end
        % Title string (to be interated so the location is by the name
        % of the tab)
            export_tabname      = sprintf('%s_%d.%d',site_tag,trans,pnt);
        
        % Export the chamber datetime objects (must import as strings to
        % allow xlswrite() to work)
        export_ChamON_times     = string(NGBN_STMB_ChamON_DateTimes(:,1,i));
        % Export the end of the flux measurement time
        export_flx_end_time     = flx_end_time;
        % Emport Chamber Data
        export_chamON_data      = NGBN_STMB_ChamON(:,2:end,i);
        % Export regression data
            % Linear (error is 95% Conf. Interval)
        export_lin_mdl_CH4      = lin_slope(:,:,i,1);
        export_lin_mdl_CO2      = lin_slope(:,:,i,2);
        export_lin_flux_CH4     = lin_flux(i,1,1);
        export_lin_flux_CH4_LB  = lin_flux(i,2,1);
        export_lin_flux_CH4_UB  = lin_flux(i,3,1);
        export_lin_flux_CO2     = lin_flux(i,1,2);
        export_lin_flux_CO2_LB  = lin_flux(i,2,2);
        export_lin_flux_CO2_UB  = lin_flux(i,3,2);
            % Exponential
        export_exp_mdl          = exp_slope(:,:,i);
        export_exp_flux         = exp_flux(i);
        % Export Keeling Plot (KP) data
        export_KP_coeffs_CH4    = KP_coeffs(:,:,i,1);
        export_KP_coeffs_CO2    = KP_coeffs(:,:,i,2);
        export_KP_coeffs_ci_CH4 = KP_coeffs_ci(:,:,i,1);
        export_KP_coeffs_ci_CO2 = KP_coeffs_ci(:,:,i,2);
        % Export last 10 isotopic composition (error is ± 1σ)
        export_L10_d13CH4       = last10_d13CH4(i,1);
        export_L10_d13CH4_LB    = last10_d13CH4(i,1) - last10_d13CH4(i,2);
        export_L10_d13CH4_UB    = last10_d13CH4(i,1) + last10_d13CH4(i,2);
        export_L10_d13CO2       = last10_d13CO2(i,1);
        export_L10_d13CO2_LB    = last10_d13CO2(i,1) - last10_d13CO2(i,2);
        export_L10_d13CO2_UB    = last10_d13CO2(i,1) + last10_d13CO2(i,2);        
        % Export Pre- and Post-backgound times (import as strings)
        export_PreB_Times       = string(NGBN_STMB_PreBack_DateTimes(:,1,i));
%         export_PostB_Times      = string(NGBN_STMB_PostBack_DateTimes(:,1,i));
        % Export Pre- and Post-background gas data
        export_PreB_CH4         = NGBN_STMB_PreBack (:,3,i);
        export_PreB_CO2         = NGBN_STMB_PreBack (:,7,i);
        export_PreB_iCH4        = NGBN_STMB_PreBack (:,5,i);
        export_PreB_iCO2        = NGBN_STMB_PreBack (:,8,i);
%         export_PostB_CH4        = NGBN_STMB_PostBack(:,3,i);
%         export_PostB_CO2        = NGBN_STMB_PostBack(:,7,i);
%         export_PostB_iCH4       = NGBN_STMB_PostBack(:,5,i);
%         export_PostB_iCO2       = NGBN_STMB_PostBack(:,8,i);
        fn                      = sprintf("%s_%s_Measurements.xlsx",    ...
                                            site_tag,ddmmmyyyy);
        XL_filename             = working_dir+fn;
        % Write the data into excel files for each variable
            % Datetimes from chamber enclosure period
        writematrix(export_ChamON_times    , XL_filename, 'Sheet', export_tabname, 'Range', 'A3');
            % Pertinent metrics from chamber enclosure period (e.g., [CH4])
        writematrix(export_chamON_data     , XL_filename, 'Sheet', export_tabname, 'Range', 'B3');
            % Linear Regression used for CH4 flux and corresponding metrics
        writematrix(export_lin_mdl_CH4     , XL_filename, 'Sheet', export_tabname, 'Range', 'M3');
            % CH4 Flux and 95% lower and upper bound
        writematrix(export_lin_flux_CH4    , XL_filename, 'Sheet', export_tabname, 'Range', 'R3');
        writematrix(export_lin_flux_CH4_LB , XL_filename, 'Sheet', export_tabname, 'Range', 'R4');
        writematrix(export_lin_flux_CH4_UB , XL_filename, 'Sheet', export_tabname, 'Range', 'R5');
            % CH4 Exponential Model
        writematrix(export_exp_mdl         , XL_filename, 'Sheet', export_tabname, 'Range', 'T3');
            % CH4 Exponential Flux
        writematrix(export_exp_flux        , XL_filename, 'Sheet', export_tabname, 'Range', 'X3');
            % Linear Regression used for CO2 flux and corresponding metrics
        writematrix(export_lin_mdl_CO2     , XL_filename, 'Sheet', export_tabname, 'Range', 'Z3');
            % CO2 Flux and 95% lower and upper bound
        writematrix(export_lin_flux_CO2    , XL_filename, 'Sheet', export_tabname, 'Range', 'AE3');
        writematrix(export_lin_flux_CO2_LB , XL_filename, 'Sheet', export_tabname, 'Range', 'AE4');
        writematrix(export_lin_flux_CO2_UB , XL_filename, 'Sheet', export_tabname, 'Range', 'AE5');
            % Keeling Plots CH4 and CO2
        writematrix(export_KP_coeffs_CH4   , XL_filename, 'Sheet', export_tabname, 'Range', 'AG3');
        writematrix(export_KP_coeffs_ci_CH4, XL_filename, 'Sheet', export_tabname, 'Range', 'AG4');
        writematrix(export_KP_coeffs_CO2   , XL_filename, 'Sheet', export_tabname, 'Range', 'AJ3');
        writematrix(export_KP_coeffs_ci_CO2, XL_filename, 'Sheet', export_tabname, 'Range', 'AJ4');
            % d13CH4 Last 10 method
        writematrix(export_L10_d13CH4      , XL_filename, 'Sheet', export_tabname, 'Range', 'AM3');
        writematrix(export_L10_d13CH4_LB   , XL_filename, 'Sheet', export_tabname, 'Range', 'AN3');
        writematrix(export_L10_d13CH4_UB   , XL_filename, 'Sheet', export_tabname, 'Range', 'AO3');
            % d13CO2 Last 10 method
        writematrix(export_L10_d13CO2      , XL_filename, 'Sheet', export_tabname, 'Range', 'AQ3');
        writematrix(export_L10_d13CO2_LB   , XL_filename, 'Sheet', export_tabname, 'Range', 'AR3');
        writematrix(export_L10_d13CO2_UB   , XL_filename, 'Sheet', export_tabname, 'Range', 'AS3');
            % Pre- and Post-Background [CH4] and [CO2]
        writematrix(export_PreB_Times      , XL_filename, 'Sheet', export_tabname, 'Range', 'AU3');
        writematrix(export_PreB_CH4        , XL_filename, 'Sheet', export_tabname, 'Range', 'AV3');
        writematrix(export_PreB_CO2        , XL_filename, 'Sheet', export_tabname, 'Range', 'AW3');
        writematrix(export_PreB_iCH4       , XL_filename, 'Sheet', export_tabname, 'Range', 'AX3');
        writematrix(export_PreB_iCO2       , XL_filename, 'Sheet', export_tabname, 'Range', 'AY3');       
%         writematrix(export_PostB_Times     , XL_filename, 'Sheet', export_tabname, 'Range', 'BA3');
%         writematrix(export_PostB_CH4       , XL_filename, 'Sheet', export_tabname, 'Range', 'BB3');
%         writematrix(export_PostB_CO2       , XL_filename, 'Sheet', export_tabname, 'Range', 'BC3');
%         writematrix(export_PostB_iCH4      , XL_filename, 'Sheet', export_tabname, 'Range', 'BD3');
%         writematrix(export_PostB_iCO2      , XL_filename, 'Sheet', export_tabname, 'Range', 'BE3');
        % Export alpha value (alpha:d13C-CO2--d13C-CH4) for Horita
        % Geothermometer
        export_alpha            = Alpha(i);
        export_alpha_Forde19    = alpha_Forde19(i);
        export_thous_ln_alpha   = thous_ln_alpha(i);
        % Write data to excel
        writematrix(export_alpha           , XL_filename, 'Sheet', export_tabname, 'Range', 'BG3');
        writematrix(export_alpha_Forde19   , XL_filename, 'Sheet', export_tabname, 'Range', 'BG4');
        writematrix(export_thous_ln_alpha  , XL_filename, 'Sheet', export_tabname, 'Range', 'BH3');
        % Message to indicate one tab has been filled in with data
        msg_str = sprintf('Tab ''%s_%d'' has been created', site_tag, i);
        disp(msg_str)
    end
        disp('All data has been exported to the worksheet')

    % Export Data to Summary Tab
    
      % Set Summary Tab  name
      export_summary_tabname = 'Summary';
      % Export Start and End Time of Enclosures
      export_SUMTAB_srt_enclosure_time =                                ...
                 reshape(string(NGBN_STMB_ChamON_DateTimes(1,1,:)), [nchams 1]);
         
      % Acquire end time index for each chamber enclosure
      j = 1;
      export_SUMTAB_end_enclosure_time = string();
      for idx = 1:nchams
          end_time_idx          = find(~isnan(NGBN_STMB_ChamON(:,1,idx)),1,  ...
                                        'last'); 
          export_SUMTAB_end_enclosure_time(j) = string(                 ...
                                      NGBN_STMB_ChamON_DateTimes(            ...
                                      end_time_idx,1,idx));
          j = j + 1;
      end
        export_SUMTAB_end_enclosure_time  =                             ...
                reshape(export_SUMTAB_end_enclosure_time,[nchams 1]);
        
        clearvars j end_time_idx

        % Arrange d13CH4 results for SUMTAB export
        export_SUMTAB_d13CH4    = reshape(KP_coeffs(:,2,:,1)   , [nchams 1]);
        export_SUMTAB_LB_d13CH4 = reshape(KP_coeffs_ci(1,2,:,1), [nchams 1]);
        export_SUMTAB_UB_d13CH4 = reshape(KP_coeffs_ci(2,2,:,1), [nchams 1]);
        export_SUMTAB_d13CO2    = reshape(KP_coeffs(:,2,:,2)   , [nchams 1]);
        export_SUMTAB_LB_d13CO2 = reshape(KP_coeffs_ci(1,2,:,2), [nchams 1]);
        export_SUMTAB_UB_d13CO2 = reshape(KP_coeffs_ci(2,2,:,2), [nchams 1]);
        
        % Write summary data into excel file
            % Total Enclosure Time
        writematrix(export_SUMTAB_srt_enclosure_time , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'H3');
        writematrix(export_SUMTAB_end_enclosure_time , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'I3');
            % Timestamps for portion used to calculate flux
        writematrix(flx_srt_time                     , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'K3');
        writematrix(flx_end_time                     , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'L3');
            % CH4 Linear Flux
        writematrix(lin_flux(:,1,1)                  , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'M3');
        writematrix(lin_flux(:,2,1)                  , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'N3');
        writematrix(lin_flux(:,3,1)                  , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'O3');
            % CH4 Keeling Plot values
        writematrix(export_SUMTAB_d13CH4             , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'P3');
        writematrix(export_SUMTAB_LB_d13CH4          , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'Q3');
        writematrix(export_SUMTAB_UB_d13CH4          , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'R3');
            % CH4 Last 10 Isotopes
        writematrix(last10_d13CH4(:,1)               , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'S3');
        writematrix(last10_d13CH4(:,1)-last10_d13CH4(:,2)               ...
                                                     , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'T3');
        writematrix(last10_d13CH4(:,1)+last10_d13CH4(:,2)               ...
                                                     , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'U3');
            % CO2 Linear Flux
        writematrix(lin_flux(:,1,2)                  , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'V3');
        writematrix(lin_flux(:,2,2)                  , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'W3');
        writematrix(lin_flux(:,3,2)                  , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'X3');
            % CO2 Keeling Plot 
        writematrix(export_SUMTAB_d13CO2             , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'Y3');
        writematrix(export_SUMTAB_LB_d13CO2          , XL_filename, 'Sheet', export_summary_tabname, 'Range', 'Z3');
        writematrix(export_SUMTAB_UB_d13CO2          , XL_filename, 'Sheet', export_summary_tabname, 'Range','AA3');
            % Last 10 CO2 Isotopes
        writematrix(last10_d13CO2(:,1)               , XL_filename, 'Sheet', export_summary_tabname, 'Range','AB3');
        writematrix(last10_d13CO2(:,1)-last10_d13CO2(:,2)               ...
                                                     , XL_filename, 'Sheet', export_summary_tabname, 'Range','AC3');
        writematrix(last10_d13CO2(:,1)+last10_d13CO2(:,2)               ...
                                                     , XL_filename, 'Sheet', export_summary_tabname, 'Range','AD3');
                    
        disp('Summary data has been written to excel sheet')
        
        disp('All data has been exported to the worksheet')
end

clearvars export_tab_name         export_ChamON_times    export_chamON_data
clearvars export_lin_flux_CH4     export_lin_mdl_CO2     export_lin_flux_CH4
clearvars export_lin_flux_CO2     export_exp_mdl         export_exp_flux
clearvars export_KP_coeffs_CH4    export_KP_coeffs_CO2   export_KP_coeffs_ci_CH4
clearvars export_KP_coeffs_ci_CO2 export_PreB_Times      
clearvars export_PreB_CH4         export_PreB_CO2        export_PreB_iCH4
clearvars export_PreB_iCO2           

clearvars export_SUMTAB_srt_time  export_SUMTAB_end_time export_SUMTAB_d13CH4
clearvars export_SUMTAB_LB_d13CH4 export_SUMTAB_UB_d13CH4 
clearvars export_SUMTAB_d13CO2    export_SUMTAB_LB_d13CO2 
clearvars export_SUMTAB_UB_d13CO2 export_summary_tabname

clearvars export_L10_d13CH4       export_L10_d13CH4_LB  export_L10_d13CH4_UB
clearvars export_L10_d13CO2       export_L10_d13CO2_LB  export_L10_d13CO2_UB
clearvars export_lin_flux_CH4_LB  export_lin_flux_CH4_UB export_lin_flux_CO2_LB
clearvars export_lin_flux_CO2_UB  export_lin_mdl_CH4

clearvars export_SUMTAB_end_enclosure_time export_SUMTAB_srt_enclosure_time

%% Export a MATFILE

% Add R-squared values from linear fits to .mat files
 Rsquared_CH4 = NaN([nchams 1]);
 Rsquared_CO2 = NaN([nchams 1]);
for i = 1:nchams
    Rsquared_CH4(i) = lin_mdl{:,i,1}.Rsquared.Ordinary;
    Rsquared_CO2(i) = lin_mdl{:,i,2}.Rsquared.Ordinary;
end
    % This will allow data to be manipulated in other scripts or elsewhere
    % in the MATLAB environment
    fn_mat   = sprintf("%s_%s_MATFILE.mat", site_tag, ddmmmyyyy);
    filename = working_dir+fn_mat;
    save(filename,                                                      ...
                  'NGBN_STMB_PreBack_DateTimes' , 'NGBN_STMB_PreBack' , ...
                  'NGBN_STMB_ChamON_DateTimes'  , 'NGBN_STMB_ChamON'  , ...
                  'lin_flux'                    , 'lin_slope'         , ...
                  'lin_mdl'                     ,                       ...
                  'Rsquared_CH4'                , 'Rsquared_CO2'      , ...
                  'last10_d13CH4'               ,                       ...
                  'last10_d13CO2'               , 'TT_NGBN_STMB'      , ...
                  '-v7.3');
                        