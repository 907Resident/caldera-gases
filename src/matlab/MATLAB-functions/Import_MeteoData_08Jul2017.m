%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Meteorological Data\Valles\DRI_WeatherData_08-13July2017.xlsx
%    Worksheet: 08Jul2017_HQ
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.


%% Import the data
[~, ~, raw0_0] = xlsread('C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Meteorological Data\Valles\DRI_WeatherData_08-13July2017.xlsx','08Jul2017_HQ','A8:A31');
[~, ~, raw0_1] = xlsread('C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Meteorological Data\Valles\DRI_WeatherData_08-13July2017.xlsx','08Jul2017_HQ','H8:H31');
[~, ~, raw0_2] = xlsread('C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Meteorological Data\Valles\DRI_WeatherData_08-13July2017.xlsx','08Jul2017_HQ','P8:P31');
[~, ~, raw0_3] = xlsread('C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Meteorological Data\Valles\DRI_WeatherData_08-13July2017.xlsx','08Jul2017_HQ','U8:U31');
[~, ~, raw0_4] = xlsread('C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Meteorological Data\Valles\DRI_WeatherData_08-13July2017.xlsx','08Jul2017_HQ','X8:X31');
raw = [raw0_0,raw0_1,raw0_2,raw0_3,raw0_4];

%% Create output variable
data = reshape([raw{:}],size(raw));

% Initiate the day of the measurements
date_of_meas = datenum('08-Jul-2017');

%% Allocate imported array to column variable names
MeteoTime   = data(:,1);
    MeteoTime   = datestr(MeteoTime + date_of_meas, 0);
    MeteoTime   = datetime(MeteoTime);
AirTemp     = data(:,2);
RelHum      = data(:,3);
BaroPress   = data(:,4);
Precip      = data(:,5);

MeteoMeas_08Jul2017   = table(MeteoTime, AirTemp, RelHum, BaroPress, Precip);  

%% Clear temporary variables
clearvars data raw raw0_0 raw0_1 raw0_2 raw0_3 raw0_4 MateoTime AirTemp RelHum BaroPress Precip