%% Project assignment 1
% The following code reads in the abalone dataset.
cdir = fileparts(mfilename('fullpath'));
% Since comma is used in the data, we need to convert these to dots
Data = fileread(fullfile(cdir,'../project/Data_abalone.csv'));
Data = strrep(Data, ',', '.');
FID = fopen('Data_abalone_dot.csv', 'w');
fwrite(FID, Data, 'char');
fclose(FID);
file_path = fullfile(cdir,'../project/Data_abalone_dot.csv');
% Table is now stored in Matlab
abalone_table = readtable(file_path);
X = table2array(abalone_table(:, 2:9));

% summary statistics for numerical attributes
S = {'Length','Diameter','Height','Whole weight','Shucked weight','Viscera weight','Shell weight','Rings'};
% summary stadistics
% the median will be studied by the percentiles/quantile. The mean is the
% quantile .5p. 
stats_names = {'Mean' 'Quantiles' 'Variance' 'Standar deviation' 'range'};

for i = 1:8
    S(i);
    m(1, i) = mean(X(1:end,i));
    q(i, :) = quantile(X(1:end,i),[0 0.25 0.5 0.75 1]);
    v(1, i) = var(X(1:end,i));
    s(1, i) = std(X(1:end,i));
    r(1, i) = range(X(1:end, i));
end

% Create a table with the stadistics summary
Stats = table(m', q, v', s', r', 'VariableNames', stats_names, 'RowNames', S');

% study the covariance and correlation between atributes
cov_matrix = cell(8, 8);
for i = 1:8
    for j=1:8
        % Calculate the covariance matrix of every attribute 
        covariance = cov(X(1:end, i), X(1:end, j));
        cov_matrix{i, j} = covariance;
        % Calculate the correlation of every attribute
        correlation(i, j) = corr(X(1:end, i), X(1:end, j));
    end
end

% What to do with the nominal?
mode(categorical(table2cell(abalone_table(:,1))))
% Sample text
