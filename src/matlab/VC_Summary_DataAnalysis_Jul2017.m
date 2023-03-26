%% Summary_Charts_VC_Jul_2017
    % Summary of Data from Jul2017 Campaign in Valles Caldera, NM
    % Five locations were measured:
        % Sulphur Springs (VCSS)
        % Alamo Canyon (VCAC)
        % Soda Dam (VCSD)
        % Redondo Peak (VCRP)
    % Metadata for these measuerments are located in each folder containing
    % raw data files
    
%% Import Data

Import_VC_Jul2017_Summary_Data

% Format timestamps for .Start_Enclosure and .End_Enclosure
VC_2017_SUMDATA.StartTimeofChamberEnclosure = datestr(                     ...
                VC_2017_SUMDATA.StartTimeofChamberEnclosure, 'HH:MM:ss');

VC_2017_SUMDATA.EndTimeofChamberClosure = datestr(                         ...
                VC_2017_SUMDATA.EndTimeofChamberClosure, 'HH:MM:ss');            

%% Group Scatter Plots
    % Group scatter plots will be useful

% Group Scatter for d13C-CO2 vs d13C-CH4
figure,
gs_isos =gscatter(VC_2017_SUMDATA.d13CH4_source,                        ...
         VC_2017_SUMDATA.d13CO2_source,                                 ...
         VC_2017_SUMDATA.Group);                              
         grid on                                                        
         xlabel('\delta^{13}C-CH_4 (‰)'), ylabel('\delta^{13}C-CO_2 (‰)')
         % Purple Diamond Marker - Sulphur Springs     = [0.85,0.70,1.00]
         set(gs_isos(4), 'Marker', 'diamond', 'Color', 'k',             ...
                         'MarkerFaceColor', [0.85,0.70,1.00])
         % Red Inverted Triangle - Alamo Canyon         = [1.00,0.40,0.40]
         set(gs_isos(1), 'Marker', 'v', 'Color', 'k',                   ...                                 
                         'MarkerFaceColor', [1.00,0.40,0.40])         
         % Teal Square - Soda Dam                       = [0.00,0.75,0.75]
         set(gs_isos(3), 'Marker', 'square', 'Color', 'k',              ...
                         'MarkerFaceColor', [0.00,0.75,0.75])         
         % Black Circle - Redondo Peak                  = [0.24,0.24,0.24]
         set(gs_isos(2), 'Marker', 'o', 'Color', 'k',                   ...
                         'MarkerFaceColor', [0.24,0.24,0.24])
        
% Group Scatter for CH4 Soil Flux and Soil Temperature
figure,
gs_CH4Flx_SoilTemp = gscatter(VC_2017_SUMDATA.SoilTemperature,          ...
                              VC_2017_SUMDATA.CH4_Flux,                 ...
                              VC_2017_SUMDATA.Group );
         grid on                                                        
         xlabel('Soil Temperature (°C)')
         ylabel('F_{CH_{4}} mg CH_4 m^{-2} hr^{-1}')
         % Purple Diamond Marker - Sulphur Springs     = [0.85,0.70,1.00]
         set(gs_CH4Flx_SoilTemp(4), 'Marker', 'diamond', 'Color', 'k',  ...
                         'MarkerFaceColor', [0.85,0.70,1.00])
         % Red Inverted Triangle - Alamo Canyon         = [1.00,0.40,0.40]
         set(gs_CH4Flx_SoilTemp(1), 'Marker', 'v', 'Color', 'k',        ...                                 
                         'MarkerFaceColor', [1.00,0.40,0.40])         
         % Teal Square - Soda Dam                       = [0.00,0.75,0.75]
         set(gs_CH4Flx_SoilTemp(3), 'Marker', 'square', 'Color', 'k',   ...
                         'MarkerFaceColor', [0.00,0.75,0.75])         
         % Black Circle - Redondo Peak                  = [0.24,0.24,0.24]
         set(gs_CH4Flx_SoilTemp(2), 'Marker', 'o', 'Color', 'k',        ...
                         'MarkerFaceColor', [0.24,0.24,0.24])
                     
                     
                     
%% Boxplots

% CH4 Fluxes
figure, 
bx_CH4_Flux = boxplot(VC_2017_SUMDATA.CH4_Flux, VC_2017_SUMDATA.Group,   ...
              'BoxStyle', 'outline','Notch', 'off');
               ylabel('F_{CH_{4}} mg CH_4 m^{-2} hr^{-1}')
               set(gca, 'YScale', 'log')
               grid on
% CO2 Fluxes
figure, 
bx_CO2_Flux = boxplot(VC_2017_SUMDATA.CO2_Flux, VC_2017_SUMDATA.Group,   ...
              'BoxStyle', 'outline','Notch', 'off');
               ylabel('F_{CO_{2}} mg CO_2 m^{-2} hr^{-1}')
               set(gca, 'YScale', 'log')
               grid on

%% Statistical Tests

