
%% Zero-Air Test

% ZAir_Aug2015_DataAnalysis
% Geoscientists: Ajayi, M.
% Location: Stevenson Center 5712 Laboratory
% Data Type: Zero air cylinder
% Metadata File: "C:\Users\moyoa\Google Drive\CompSci\MS Thesis\Data\Collected_Data\Nashville\Picarro\August2015\25Aug2015\Picarro_Metadata_25Aug2015.txt"

%% Extra Functions Used in this Script
    % The following functions that were made by Ajayi, M. Their
    % documentation should be located in the "Functions" folder under
    % <"R:/Ayers/FrackingTN/Data/Matlab/Data Analysis">
    
        % lin_regress()
        % KeelingCurve()
        % StatsPlot_lite()
        % CH4_vs_Time()
        % CO2_vs_Time()
        % hp_filt()
        % fullscreen()

%% Gather the ZAir_Aug2015 Data from the Forerunner Matrix 'PD_mtrx'

% Time Data (fraction of a day)
ZAir_Aug2015_X1         = PD_mtrx(:, 1);        % unitless
% Time Data (seconds relative to the start of the measuring period)
     % Change the serial time(fraction of a day) into seconds (the second
     % of the day, anywhere between 1 and 86400
ZAir_Aug2015_X1_sec     = 24 * 3600 .* ZAir_Aug2015_X1;   % seconds
    % Change the seconds to make them relative to t_0 (from 0 to 2400, or
    % 40 min)
ZAir_Aug2015_X1_sec     = ZAir_Aug2015_X1_sec - ...
                              min(ZAir_Aug2015_X1_sec);

% Important Metrics
ZAir_Aug2015_HP_12CH4_dry       = PD_mtrx(:, 18);        % ppm
ZAir_Aug2015_HP_Delta_iCH4_Raw  = PD_mtrx(:, 22);        % ppm
ZAir_Aug2015_Delta_Raw_iCO2     = PD_mtrx(:, 38);        % ppm
ZAir_Aug2015_CO2_dry            = PD_mtrx(:, 35);        % ppm

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
    % For the ZAir_Aug2015, ZAir_Aug2015, and BHCS sites, the data presented groups of four
    % thus, only every fourth row is retainted
RN = zeros([length(ZAir_Aug2015_HP_Delta_iCH4_Raw), 1]);
    % Assigns '1' to every fourth row in the vector RN
for i = 1:4:length(RN)
RN(i) = 1;
end
    % Create a new matrix that will contain the desired measurements and
    % the keep/delete indicator by horizontally concatenating RN to the
    % measurement vectors
RN_CH4 = [RN, ZAir_Aug2015_X1, ZAir_Aug2015_HP_12CH4_dry, ZAir_Aug2015_HP_Delta_iCH4_Raw, ...
          ZAir_Aug2015_Delta_Raw_iCO2, ZAir_Aug2015_CO2_dry];
    % Change the soon-to-be deleted rows to NaNs
for i = 1:length(RN)
    if RN_CH4(i,1) == 0
       RN_CH4(i,:) = NaN;
    end
    %*%* This line should be probably be deleted %*%*
    ZAir_Aug2015_X1(i) = RN_CH4(i,2);
end
    % Delete the rows that contain NaNs (as determined by the indicator
    % vector)
RN_CH4(any(isnan(RN_CH4),2),:) = [];

% Reassign the variables so that they represent the filtered information 
    ZAir_Aug2015_X1                 = RN_CH4(:,2);
    ZAir_Aug2015_X1_sec             = 24 * 3600 .* ZAir_Aug2015_X1;
    ZAir_Aug2015_X1_sec             = ZAir_Aug2015_X1_sec - ...
                                          min(ZAir_Aug2015_X1_sec);
    ZAir_Aug2015_HP_12CH4_dry       = RN_CH4(:,3);
    ZAir_Aug2015_HP_Delta_iCH4_Raw  = RN_CH4(:,4);
    ZAir_Aug2015_Delta_Raw_iCO2     = RN_CH4(:,5);
    ZAir_Aug2015_CO2_dry            = RN_CH4(:,6);

% Visualize the Data
%------Figure 01-------%
StatsPlot_lite(ZAir_Aug2015_X1, ZAir_Aug2015_HP_12CH4_dry, 13)
    title('Zero Air Test Aug2015 [CH_4] vs. Timestamp (ALL)')
    ylabel('[CH_4] (ppm)')
    legend('[CH_4]', '± 1 std', 'Mean', 'Max', 'Min', ...
           'Location', 'Northwest')
       
%------Figure 02-------%
StatsPlot_lite(ZAir_Aug2015_X1, ZAir_Aug2015_CO2_dry, 13)
    title('Zero Air Test Aug2015 [CO_2] vs. Timestamp (ALL)')
    ylabel('[CO_2] (ppm)')
    legend('[CO_2]', '± 1 std', 'Mean', 'Max', 'Min', ...
           'Location', 'Northwest')
       
%------Figure 03-------%
StatsPlot_lite(ZAir_Aug2015_X1(977:1445), ZAir_Aug2015_CO2_dry(977:1445), 13)
    title('Zero Air Test Aug2015 [CO_2] vs. Timestamp (Zero Air ONLY)')
    ylabel('[CO_2] (ppm)')
    legend('[CO_2]', '± 1 std', 'Mean', 'Max', 'Min', ...
           'Location', 'Northwest')
        