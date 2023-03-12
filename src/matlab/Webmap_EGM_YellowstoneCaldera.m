%% Webmap Valles Caldera

%% Import Data

Import_EGM_YNP_2018

%% Create Webmap for Valles Caldera
    % Use the USGS Imagery Webmap 
    wm_YC = webmap('USGS Imagery');
    % Set up axis limits to narrow the scope of the webmap
        % These hard coded coordinates are chosen upon visual inspection
    latlim = [45.102589, 44.124510]; lonlim = [-111.159435, -109.801318];
    wmlimits(wm_YC, latlim, lonlim)
    % Zoom the webmap toward the caldera
    wmzoom(wm_YC, 12)
    
 %% Create Geopoint Struct
    % The geopoint struct will store the information needed to eventually
    % create the webmap markers and the kml file

% Create an empty geopoint structure
YC_EGM_GeoPts = geopoint();

% Define the geometry of the geopoint struct
YC_EGM_GeoPts.Geometry = deal('Point');

% Extract the Latitude & Longitude and enter it into the geopoint struct
YC_EGM_GeoPts.Latitude  = EGM_YNP_SUMDATA_2018.EGM_Lat;
YC_EGM_GeoPts.Longitude = EGM_YNP_SUMDATA_2018.EGM_Lon;

% Add the names to each of the points
YC_EGM_GeoPts.Name = EGM_YNP_SUMDATA_2018.Location;

% Create a description variable (w/ HTML tags) that will include the
% relevant measurement values
description = string();

% Loop through the various points to define the description variable for
% each point in the struct
for i = 1:height(EGM_YNP_SUMDATA_2018)
description(i) = sprintf("<ul>  <li>CO2 Flux = %3.3e g m-2 hr-1</li> <li>log10[CO2 Flux] = %3.3e g m-2 hr-1</li> <li>Dist from Picarro = %2.2f m</li> <li>Location = %2.2f per mil</li> </ul>", ...
                      EGM_YNP_SUMDATA_2018.fCO2(i),                     ...
                      EGM_YNP_SUMDATA_2018.log10fCO2(i),                ...
                      EGM_YNP_SUMDATA_2018.Distance_from_Picarro_m(i),  ...
                      EGM_YNP_SUMDATA_2018.Site(i));
end

YC_EGM_GeoPts.Description = description;

clearvars description
%% Plot Points in Webmap
    % Use the webmap handle created earlier and load in the geopint struct
    % as points into the map
    
wmmarker(wm_YC, YC_EGM_GeoPts, 'FeatureName', YC_EGM_GeoPts.Name)

%% Export Geopoint Struct to KML file

% Define the filename for the kml file
kml_filename = "C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Picarro\Methodology\Maps & Geography\New Mexico\YC_EGM_Jul2017.kml";

kmlwritepoint(kml_filename, YC_EGM_GeoPts.Latitude, YC_EGM_GeoPts.Longitude,            ...
             'Description', YC_EGM_GeoPts.Description, 'Name', YC_EGM_GeoPts.Name,      ...
             'Icon', 'http://maps.google.com/mapfiles/kml/shapes/placemark_square.png', ...
             'Color', [0.4844 0.0391 0.0078],                                           ...
             'IconScale', 0.75)
