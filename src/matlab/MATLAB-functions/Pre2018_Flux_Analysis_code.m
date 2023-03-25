%% Pre2018 Flux Analysis

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
       CH4_lin_mdl{i}           = fitlm(VCAC_11Jul2017_ISM_DATA(:,1,i),           ...
                                        VCAC_11Jul2017_ISM_DATA(:,8,i));
       CH4_lin_slope(:,:,i)     = table2array(CH4_lin_mdl{i}.Coefficients);
       % Quantify the flux for each measurement (mg m^-2 hr^-1)
       CH4_lin_flux(i)          = CH4_lin_slope(2,1,i) .* H .* hour;
    end

        %----Carbon Dioxide----%
    % Pre-allocate vector for fluxes (per measurement cycle)
    CO2_lin_mdl     = cell ([1 nchams]);
    CO2_lin_slope   = zeros([2 4 nchams]);
    CO2_lin_flux    = zeros([1 nchams]);
    % Use for-loop to run a linear regression for each of the chamber
    % measuerments
    for i = 1:nchams
       CO2_lin_mdl{i}           = fitlm(VCAC_11Jul2017_ISM_DATA(:,1,i), ...
                                        VCAC_11Jul2017_ISM_DATA(:,9,i));
       CO2_lin_slope(:,:,i)     = table2array(CO2_lin_mdl{i}.Coefficients);
       % Quantify the flux for each measurement (mg m^-2 hr^-1)
       CO2_lin_flux(i)          = CO2_lin_slope(2,1,i) .* H .* hour;
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
                        prepareCurveData(VCAC_11Jul2017_ISM_DATA(:,1,i), ...
                                         VCAC_11Jul2017_ISM_DATA(:,8,i));
x_expData{i}            = x_exp_prepData;
y_expData{i}            = y_exp_prepData;
CH4_exp_mdl{i}          = fit(x_expData{i},y_expData{i},CH4_exp_fxn,              ...
                             'problem', VCAC_11Jul2017_ISM_DATA(1,8,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
CH4_exp_slope(1,1,i)  = CH4_exp_mdl{i}.psi;
CH4_exp_slope(1,2,i)  = CH4_exp_mdl{i}.kappa;
CH4_exp_slope(1,3,i)  = CH4_exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
CH4_exp_flux(i)       =((CH4_exp_mdl{i}.psi - CH4_exp_mdl{i}.c_0) .* CH4_exp_mdl{i}.kappa) .* H .* hour;
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
                        prepareCurveData(VCAC_11Jul2017_ISM_DATA(:,1,i), ...
                                         VCAC_11Jul2017_ISM_DATA(:,9,i));
x_expData{i}            = x_exp_prepData;
y_expData{i}            = y_exp_prepData;
CO2_exp_mdl{i}          = fit(x_expData{i},y_expData{i},CO2_exp_fxn,              ...
                             'problem', VCAC_11Jul2017_ISM_DATA(1,9,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
CO2_exp_slope(1,1,i)  = CO2_exp_mdl{i}.psi;
CO2_exp_slope(1,2,i)  = CO2_exp_mdl{i}.kappa;
CO2_exp_slope(1,3,i)  = CO2_exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
CO2_exp_flux(i)       =((CO2_exp_mdl{i}.psi - CO2_exp_mdl{i}.c_0) .* CO2_exp_mdl{i}.kappa) .* H .* hour;
end
clearvars x_expData y_expData x_exp_prepData y_exp_prepData
