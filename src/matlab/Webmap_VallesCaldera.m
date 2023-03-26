%% Webmap Valles Caldera

%% Import Data

Import_VC_2017_Summary_Data

%% Create Webmap for Valles Caldera
    % Use the USGS Imagery Webmap 
    wm_VC = webmap('USGS Imagery');
    % Set up axis limits to narrow the scope of the webmap
        % These hard coded coordinates are chosen upon visual inspection
    latlim = [35.7857, 36.021696]; lonlim = [-106.8 -106.4];
    wmlimits(wm_VC, latlim, lonlim)
    % Zoom the webmap toward the caldera
    wmzoom(wm_VC, 12)
    
 %% Create Geopoint Struct
    % The geopoint struct will store the information needed to eventually
    % create the webmap markers and the kml file

% Create an empty geopoint structure
VC_GeoPts = geopoint();
    
% Define the geometry of the geopoint struct
VC_GeoPts.Geometry = deal('Point');

% Extract the Latitude & Longitude and enter it into the geopoint struct
VC_GeoPts.Latitude  = VC_2017_Summary.Latitude;
VC_GeoPts.Longitude = VC_2017_Summary.Longitude;

% Add the names to each of the points
VC_GeoPts.Name = VC_2017_Summary.Location;

% Create a description variable (w/ HTML tags) that will include the
% relevant measurement values
description = string();

% Loop through the various points to define the description variable for
% each point in the struct
for i = 1:height(VC_2017_Summary)
description(i) = sprintf("<ul>  <li>CH4 Flux = %3.3e mg m-2 hr-1</li> <li>CO2 Flux = %3.3e mg m-2 hr-1</li> <li>d13C-CH4 = %2.2f per mil</li> <li>d13C-CO2 = %2.2f per mil</li> <li>Group: %s</li> </ul>", ...
                      VC_2017_Summary.CH4Flux(i),                          ...
                      VC_2017_Summary.CO2Flux(i),                          ...
                      VC_2017_Summary.d13CH4_source(i),                    ...
                      VC_2017_Summary.d13CO2_source(i),                    ...
                      VC_2017_Summary.Group(i));
end

VC_GeoPts.Description = description;

clearvars description

%% Plot Points in Webmap
    % Use the webmap handle created earlier and load in the geopint struct
    % as points into the map
    
wmmarker(wm_VC, VC_GeoPts, 'FeatureName', VC_GeoPts.Name)

%% Export Geopoint Struct to KML file

% Define the filename for the kml file
kml_filename = "C:\Users\moyoa\Google Drive\CompSci\PhD Dissertation\Data Analysis\Picarro\Methodology\Maps & Geography\New Mexico/VC_Jul2017.kml";

kmlwritepoint(kml_filename, VC_GeoPts.Latitude, VC_GeoPts.Longitude,       ...
             'Description', VC_GeoPts.Description, 'Name', VC_GeoPts.Name, ...
             'Icon', 'http://maps.google.com/mapfiles/kml/shapes/target.png', ...
             'Color', [0.4844 0.0391 0.0078])
