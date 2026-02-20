function stats_results_tbl = desc_stats(data)
%desc_stats Quantifies summary statistics (e.g., mean, median, standard
%deviation, skewness, and kurtosis of a column or pages of columns of
%numeric data in matrix form. 
%   Detailed explanation goes here

% Create function handle to conduct analyses
summary_stats    = @(x)[nanmean(x) nanmedian(x) nanstd(x) skewness(x)                    ...
                       kurtosis(x)];
% Apply the function handle to the data
stats_results    = summary_stats(data); 

% Convert results into a table
stats_results_tbl = array2table(stats_results,                                           ...
                                 'VariableNames', ["Mean", "Median", "Std. Deviation",   ...
                                                   "Skewness", "Kurtosis"]);
end

