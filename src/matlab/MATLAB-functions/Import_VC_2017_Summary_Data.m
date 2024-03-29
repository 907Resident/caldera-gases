%% Import VC 2017 Summary data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Picarro\Valles\Summary Data\VallesCaldera_Jul2017_DataSummary.xlsx
%    Worksheet: Data
%
% Auto-generated by MATLAB on 30-Aug-2019 10:06:52

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 26);

% Specify sheet and range
opts.Sheet = "Data";
opts.DataRange = "A3:Z27";

% Specify column names and types
opts.VariableNames = ["Group", "Location", "Soil_Classification", "Latitude", "Longitude", "Date_of_Measurement", "Start_Time_of_Chamber_Enclosure", "End_Time_of_Chamber_Enclosure", "Duration_of_Chamber_Enclosure", "CH4_Flux", "Long_TermCH4_Flux", "d13CH4_source", "LowerBound_d13CH4_source", "UpperBound_d13CH4_source", "CO2_Flux", "Long_Term_CO2_Flux", "d13CO2_source", "LowerBound_d13CO2_source", "UpperBound_d13CO2_source", "Horita_Geothermometer_Temperature_at_Formation", "Ambient_Temperature", "Barometric_Pressure", "Soil_Moisture", "Soil_Conductivity", "Soil_Tempeature_at_Surface", "Soil_Tempeature_at_6_in_depth"];
opts.SelectedVariableNames = ["Group", "Location", "Soil_Classification", "Latitude", "Longitude", "Date_of_Measurement", "Start_Time_of_Chamber_Enclosure", "End_Time_of_Chamber_Enclosure", "Duration_of_Chamber_Enclosure", "CH4_Flux", "Long_TermCH4_Flux", "d13CH4_source", "LowerBound_d13CH4_source", "UpperBound_d13CH4_source", "CO2_Flux", "Long_Term_CO2_Flux", "d13CO2_source", "LowerBound_d13CO2_source", "UpperBound_d13CO2_source", "Horita_Geothermometer_Temperature_at_Formation", "Ambient_Temperature", "Barometric_Pressure", "Soil_Moisture", "Soil_Conductivity", "Soil_Tempeature_at_Surface", "Soil_Tempeature_at_6_in_depth"];
opts.VariableTypes = ["categorical", "string", "string", "double", "double", "datetime", "datetime", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, 6, "InputFormat", "");
opts = setvaropts(opts, 7, "InputFormat", "");
opts = setvaropts(opts, 8, "InputFormat", "");
opts = setvaropts(opts, [2, 3, 26], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 26], "EmptyFieldRule", "auto");

% Import the data
VC_SUMDATA_2017 = readtable("C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Picarro\Valles\Summary Data\VallesCaldera_Jul2017_DataSummary.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts