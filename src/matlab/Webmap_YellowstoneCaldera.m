%% Webmap Valles Caldera

%% Import Data

Import_YC_2018_Summary_Data

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
YC_GeoPts = geopoint();

% Define the geometry of the geopoint struct
YC_GeoPts.Geometry = deal('Point');

% Extract the Latitude & Longitude and enter it into the geopoint struct
YC_GeoPts.Latitude  = YNP_SUMDATA_2018.Latitude;
YC_GeoPts.Longitude = YNP_SUMDATA_2018.Longitude;

% Add the names to each of the points
YC_GeoPts.Name = YNP_SUMDATA_2018.Location;

% Create a description variable (w/ HTML tags) that will include the
% relevant measurement values
description = string();

% Loop through the various points to define the description variable for
% each point in the struct
for i = 1:height(YNP_SUMDATA_2018)
description(i) = sprintf("<ul>  <li>CH4 Flux = %3.3e mg m-2 hr-1</li> <li>CO2 Flux = %3.3e mg m-2 hr-1</li> <li>d13C-CH4 = %2.2f per mil</li> <li>d13C-CO2 = %2.2f per mil</li> <li>Group: %s</li> </ul>", ...
                      YNP_SUMDATA_2018.CH4_Flux(i),                          ...
                      YNP_SUMDATA_2018.CO2_Flux(i),                          ...
                      YNP_SUMDATA_2018.d13CH4_source(i),                    ...
                      YNP_SUMDATA_2018.d13CO2_source(i),                    ...
                      YNP_SUMDATA_2018.Group(i));
end

YC_GeoPts.Description = description;

clearvars description
%% Plot Points in Webmap
    % Use the webmap handle created earlier and load in the geopint struct
    % as points into the map
    
wmmarker(wm_YC, YC_GeoPts, 'FeatureName', YC_GeoPts.Name)

%% Export Geopoint Struct to KML file

% Define the filename for the kml file
kml_filename = "C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Picarro\Methodology\Maps & Geography\New Mexico/YC_Jul2017.kml";

kmlwritepoint(kml_filename, YC_GeoPts.Latitude, YC_GeoPts.Longitude,       ...
             'Description', YC_GeoPts.Description, 'Name', YC_GeoPts.Name, ...
             'Icon', 'http://maps.google.com/mapfiles/kml/shapes/target.png', ...
             'Color', [0.4844 0.0391 0.0078])
