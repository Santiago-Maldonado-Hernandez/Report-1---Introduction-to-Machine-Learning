%% Project assignment 1
clear all
clc
%% Initial Commands
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
Sr = {'Length','Diameter','Height','Whole weight','Shucked weight','Viscera weight','Shell weight','Rings'};
% summary stadistics
% the median will be studied by the percentiles/quantile. The median is the
% quantile .5p. 
stats_names = {'Mean' 'Quantiles' 'Variance' 'Standar deviation' 'range'};

% iniziate vectors
m = zeros(1, 8);
q = zeros(8, 5);
v = zeros(1, 8);
s = zeros(1, 8);
r = zeros(1, 8);

for i = 1:8
    m(1, i) = mean(X(1:end,i));
    q(i, :) = quantile(X(1:end,i),[0 0.25 0.5 0.75 1]);
    v(1, i) = var(X(1:end,i));
    s(1, i) = std(X(1:end,i));
    r(1, i) = range(X(1:end, i));
end

% Create a table with the stadistics summary
Stats = table(m', q, v', s', r', 'VariableNames', stats_names, 'RowNames', Sr');

% study the covariance and correlation between atributes
cov_matrix = cell(8, 8);
correlation = zeros(8, 8);
for i = 1:8
    for j=1:8
        % Calculate the covariance matrix of every attribute 
        covariance = cov(X(1:end, i), X(1:end, j));
        cov_matrix{i, j} = covariance;
        % Calculate the correlation of every attribute
        correlation(i, j) = corr(X(1:end, i), X(1:end, j));
    end
end


%% Visualization

% Subtract the mean from the data
Y = bsxfun(@minus, X, mean(X));

% Obtain the PCA solution by calculate the SVD of Y
[U, S, V] = svd(Y);

% Compute variance explained
rho = diag(S).^2./sum(diag(S).^2);
threshold = 0.95;

% Plot variance explained
mfig('NanoNose: Var. explained'); clf;
hold on
plot(rho, 'x-');
plot(cumsum(rho), 'o-');
plot([0,length(rho)], [threshold, threshold], 'k--');
legend({'Individual','Cumulative','Threshold'}, ...
        'Location','best');
ylim([0, 1]);
xlim([1, length(rho)]);
grid minor
xlabel('Principal component');
ylabel('Variance explained value');
title('Variance explained by principal components');

% boxplots of the dimensions and different weights
figure(2)
boxplot(X(1:end,1:3),{'Length','Diameter','Height'})
xlabel('Attribute')
ylabel('Size (mm)')
figure(3)
boxplot(X(1:end,4:7),{'Whole weight','Shucked weight','Viscera weight','Shell weight'})
xlabel('Attribute')
ylabel('Weight (g)')

% Histogram plots
for i=1:8
   figure(i+3)
   H(i) = histogram(X(:, i));
   xlabel(Sr(i));
   title(['Histogram of the attribute ' Sr(i)])
end

% Plots distinguising the sex (M) or (F) or (I)
sex = table2array(abalone_table(:, 1));
figure(12)
for i=1:4177
   if(strcmp(sex{i, 1}, 'M'))
       plot(X(i, 1), X(i,4), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
       hold on
   elseif (strcmp(sex{i, 1}, 'F'))
       plot(X(i, 1), X(i,4), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
       hold on
   else 
       plot(X(i, 1), X(i,4), 'o', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y')
       hold on
   end
end
legend('Male', 'Female', 'Infant')
xlabel('Length (mm)')
ylabel('Whole weight')
title('Length vs Whole wieght taking into account the sex')
hold off

%% 
% What to do with the nominal?
mode(categorical(table2cell(abalone_table(:,1))))
% Sample text
