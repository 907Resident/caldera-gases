function [TT_PicData,PD_mtrx_nm,PD_mtrx_dt, dte] =                         ...
          Merge_PicData_Fxn(directory, TZ_Local)
%MERGEPICDATA_FXN Concatenate all .dat files from a day in the field with
%the Picarro and eosAC
 % Merge_PicData
    % This script assists in the importation and concantenation of multiple
    % .dat files containing Picarro data.
 % INPUTS:
 % - directory: 
    % It is preferable to set up the directory as:
    % ...\subfolder1\subfolder2\*.dat
 % - TZ_Local:
    % Enter the time zone of where the measurements took place. 
    % EST = America/New_York, CST = America/Chicago, MST = America/Denver, 
    % PST = America/Los_Angeles
 % OUTPUTS:
 % - TT_PicData:
    % Timeseries object that contains the concatenated .dat files, indexed
    % by time
 % - PD_mrtx_nm:
    % A matrix of double precision values that contains all of the numeric 
    % data collected in the .dat files
 % - PD_mrtx_dt:
    % A matrix of double precision values that contains all of the
    % timestamped data collected in the .dat files
%
%
% 
 % Instructions for use
   % Step 01: Call for the initial variables
        % The user informs the algorithm on how many workbooks that will be
        % imported.  **Also, the user will need to record and enter a
        % vector that indicates  the number of rows that include data in
        % each workbook.** Please read the notes on formating the Excel
        % Workbooks below to ensure the importing process works smoothly. 
%
%     
   % Step 02: Concatenate data from the excel files into one master matrix 
        % Using a 'for-loop', the algorithm  uses the filename (the file
        % path of the workbook) and imports the data into the workspace
        % using another script entitled, 'importPicDat' (for more
        % information on this script type >>> edit 'importPicDat' into
        % the Command Window).  The 'for-loop' will continue to import the
        % Picarro data workbook by workbook, and each workbook will be
        % concatenated vertically so that final product is *one* large
        % matrix entitled, 'PD_mtrx.'  Each workbook is cataloged
        % individually in a 1 x Nfiles cellular array entitled, 'PicData.'
 %  
 %------------------------------------------------------------------------%
 %
 % Created by M. Ajayi (Sep 2019)
 
%% Step 01: Call for the initial variables
    % Set up the directory to glean all of the .dat files in the directory.
D           = dir(directory);
Nfiles      = length(D);
PicData     = cell(1, Nfiles);

%% Step 02: Concatenate data from the excel files into one master matrix
    % This section goes through each of the workbook files and pulls the
    % data from each file and concatenates it into one master matrix called
    % 'PD_mtrx.'

% Prepare blank matrices for appending data in the for-loop
PD_mtrx_dt = [];
PD_mtrx_nm = [];
    
disp('Starting concatenation ...')

for i = 1:Nfiles
% Start the timer to show the length of time needed to process .dat
% file
    tic
% Load in the file name from the directory
    filename = D(i).name;

% Store all variables from .dat files in a table (use importPicDat)
    tableout            = importPicDat(filename);            
 % Separate tables by variable type

    % Gather the dates from the timestamps
    dte                 = tableout.DATE;
    dtm                 = datetime(tableout.TIME,'InputFormat',         ...
                                  'HH:mm:ss.SSS');
    dte.Format          = 'dd-MMM-yyyy HH:mm:ss.SSS';
    dtm.Format          = 'dd-MMM-yyyy HH:mm:ss.SSS';
    
    % Collect the time from the timestamps
    DateTime            = dte + timeofday(dtm);
    DateTime.TimeZone   ='UTC';
    DateTime.TimeZone   = TZ_Local;
    % Insert the datetime into the dummy table
    tableout            = addvars(tableout, DateTime, 'After', 'TIME');
    % Format the time portion of the date timestamp
    DateTime.Format     ='HH:mm:ss.SSS';
    TIME_TZ_Local       = DateTime;
   
    tableout            = addvars(tableout, TIME_TZ_Local, 'After',     ...
                                 'DateTime');
 
% Convert the table to an array
    % Dates and Times
    tableout_datestimes = tableout.DateTime;
    dmy_datestimes      = tableout_datestimes;
    
    % Numeric
    vartype_nm          = vartype('numeric');
    tableout_numeric    = tableout(:,vartype_nm);
    tableout_numeric    = table2array(tableout_numeric);
    dmy_numeric         = tableout_numeric;
     
 
    PD_mtrx_dt          = [PD_mtrx_dt; dmy_datestimes]; 
    PD_mtrx_nm          = [PD_mtrx_nm; dmy_numeric];  
                 
    PicData{i}          = PD_mtrx_nm;

    PicData{i}          = {PD_mtrx_dt, PicData{i}};
    
    toc
    
    % Progress Bar
    if i == 1
        wb = waitbar(0,'Beginning concatenation of .dat files');
             pause(0.75)
    elseif i == ceil(Nfiles * 0.15)
             waitbar(0.15,wb, '15 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.25)
             waitbar(0.25,wb, '25 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.35)
             waitbar(0.35,wb, '35 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.50)
             waitbar(0.50,wb, '50 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.65)
             waitbar(0.65,wb, '65 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.75)
             waitbar(0.75, wb, '75 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.85)
             waitbar(0.85, wb, '85 percent of concatenations complete')
             pause(0.75)
    elseif i == ceil(Nfiles * 0.95)
            waitbar(0.95, wb,  '95 percent of concatenations complete')
            pause(0.75)
    end
    
    % Close progress bar
    if i == Nfiles-1
        waitbar(0.99, wb, 'Last iteration')
    elseif i == Nfiles
        close(wb)
    end
    
end

disp('... Concatenation Complete')

% Create a Timetable object to include the pertinent variables
    % Set up variables names
    VariableNames = {'HP_CH4_dry', 'HR_CH4_dry', 'HP_d13_CH4',          ...
                     'HR_d13_CH4', 'CO2_dry', 'd13_CO2', 'Alarm_Status'};
                 
    % Create timetable
    TT_PicData = timetable(PD_mtrx_dt(:,1), PD_mtrx_nm(:,16),           ...
               PD_mtrx_nm(:,17), PD_mtrx_nm(:,20), PD_mtrx_nm(:,21),    ...
               PD_mtrx_nm(:,34), PD_mtrx_nm(:,36), PD_mtrx_nm(:,31),    ...
               'VariableNames', VariableNames);
           
disp('Timetable object created')    

clearvars dmy_numeric dmy_datestimes vartype_dt vartype_nm tableout
clearvars tableout_datestimes
clearvars dmy_numeric dmy_datestimes vartype_dt vartype_nm tableout
clearvars tableout_numeric tableout_datestimes

end