close all; clear variables; clc;

% Directory setup for saving outputs
if ~exist('results', 'dir')
    mkdir('results');
end

% Import and prepare dataset
data = importdata("superconduct.csv");
[m, n] = size(data);

preproc = 1;
[Dtrn, Dval, Dchk] = split_scale(data, preproc);

% Define kFold, Radius and Features
kFold = 5;
features = [5, 10, 15, 20];
radius = [0.3, 0.5, 0.7, 0.9];

n1 = length(features);
n2 = length(radius);

meanError = zeros(n1, n2);
numOfRules = zeros(1, n2);

rules = zeros(kFold, 1);
errors = zeros(kFold, 1);

[indexList, weights] = relieff(Dtrn(:, 1:end - 1), Dtrn(:, end), 10);

% Grid Search
for i = 1: 1: n1
    for j = 1: 1: n2
        cvpart = cvpartition(size(Dtrn, 1), 'KFold', kFold);
        for k = 1: 1: kFold

            % chosen is the new train dataset after relieff is applied on the initial train dataset 
            chosen = [Dtrn(:, indexList(1:features(i))) Dtrn(:, end)];

            % find the new indexes for train and validation dataset 
            validationIndex = find(test(cvpart, k) == 1);
            trainingIndex = find(training(cvpart, k) == 1);

            % Get the new train and validation dataset
            DtrnNew = chosen(trainingIndex, :);
            DvalNew = chosen(validationIndex, :);

            % Create Fuzzy Inference System (FIS)
            options = genfisOptions('SubtractiveClustering');
            options.ClusterInfluenceRange = radius(j);
            tsk = genfis(DtrnNew(:, 1:end - 1), DtrnNew(:, end), options);

            % Train Fuzzy Inference System (FIS)
            options2 = anfisOptions;
            options2.InitialFIS = tsk;
            options2.ValidationData = DvalNew;
            options2.EpochNumber = 100;
            options2.OptimizationMethod = 1;

            [fis, trainError, stepSize, valFIS, valError] = anfis(DtrnNew, options2);

            rules(k) = size(showrule(valFIS), 1);
            errors(k) = abs(mean(valError));
        end
        meanError(i, j) = mean(errors);
        numOfRules(j) = round(mean(rules));
    end
end

% Plots Mean Absolute Error - Number of Rules and Mean Absolute Error - Features
figure;
hold on
title('Error - Number Of Features', 'Interpreter', 'latex');
xlabel('Number Of Features', 'Interpreter', 'latex');
ylabel('Mean Error', 'Interpreter', 'latex');
plot(features, meanError(:, 1), 'o', 'LineWidth', 1, 'Color', 'red');
plot(features, meanError(:, 2), 'o', 'LineWidth', 1, 'Color', 'blue');
plot(features, meanError(:, 3), 'o', 'LineWidth', 1, 'Color', 'green');
plot(features, meanError(:, 4), 'o', 'LineWidth', 1, 'Color', 'black');
labels = num2str([numOfRules(1), numOfRules(2), numOfRules(3), numOfRules(4)].', '%d Rules');
legend(labels);
xlim([0 30]);
grid on
saveas(gcf, 'results/Error_vs_Features.png');

figure;
hold on
title('Error - Number Of Rules', 'Interpreter', 'latex');
xlabel('Number Of Rules', 'Interpreter', 'latex');
ylabel('Mean Error', 'Interpreter', 'latex');
plot(numOfRules, meanError(1, :), 'o', 'LineWidth', 1, 'Color', 'red');
plot(numOfRules, meanError(2, :), 'o', 'LineWidth', 1, 'Color', 'blue');
plot(numOfRules, meanError(3, :), 'o', 'LineWidth', 1, 'Color', 'green');
plot(numOfRules, meanError(4, :), 'o', 'LineWidth', 1, 'Color', 'black');
legend('5 features', '10 features', '15 features', '20 features');
grid on
xlim([0 15]);
saveas(gcf, 'results/Error_vs_Rules.png');

% Get Best Model
[bestIndex1, bestIndex2] = find(meanError == min(min(meanError)));
bestFeatures = features(bestIndex1);
bestRadius = radius(bestIndex2);

DtrnBest = [Dtrn(:, indexList(1:bestFeatures)) Dtrn(:, end)];
DvalBest = [Dval(:, indexList(1:bestFeatures)) Dval(:, end)];
DchkBest = [Dchk(:, indexList(1:bestFeatures)) Dchk(:, end)];

% Create Best Fuzzy Inference System (FIS)
optionsBest = genfisOptions('SubtractiveClustering');
optionsBest.ClusterInfluenceRange = bestRadius;
tskBest = genfis(DtrnBest(:, 1:end - 1), DtrnBest(:, end), optionsBest);

% Train Best Fuzzy Inference System (FIS)
optionsBest2 = anfisOptions;
optionsBest2.InitialFIS = tskBest;
optionsBest2.ValidationData = DvalBest;
optionsBest2.EpochNumber = 100;
optionsBest2.OptimizationMethod = 1;

[fisBest, trainErrorBest, stepSizeBest, valFISBest, valErrorBest] = anfis(DtrnBest, optionsBest2);

% Evaluate Best Fuzzy Inference System (FIS)
outputBest = evalfis(DchkBest(:, 1:end - 1), valFISBest);

% Plot Prediction Values and Real Values
figure;
hold on;
title('Prediction Values', 'Interpreter', 'latex');
plot(1: length(outputBest), outputBest, 'bo');
saveas(gcf, 'results/Predictions.png');

figure;
hold on;
title('Real Values', 'Interpreter', 'latex');
plot(1: length(outputBest), Dchk(:, end), 'r*');
saveas(gcf, 'results/RealValues.png');

% Plot Learning Curve
figure;
plot([trainErrorBest valErrorBest]);
legend('Training Error', 'Validation Error');
title("Learning Curve", 'interpreter','latex');
saveas(gcf, 'results/LearningCurve.png');

% Plot Membership Functions before and after training
figure;
subplot(2, 1, 1);
plotmf(tskBest, 'input', 1);
title("feature 1 before training", 'interpreter','latex');
subplot(2, 1, 2);
plotmf(valFISBest, 'input', 1);
title("feature 1 after tarining",'interpreter','latex');
saveas(gcf, 'results/MF_Input1.png');

figure;
subplot(2, 1, 1);
plotmf(tskBest, 'input', 2);
title("feature 2 before training", 'interpreter','latex');
subplot(2, 1, 2);
plotmf(valFISBest, 'input', 2);
title("feature 2 after tarining",'interpreter','latex');
saveas(gcf, 'results/MF_Input2.png');

% Calculate RMSE, NMSE, NDEI, R^2
rmse = sqrt(mse(Dchk(:, end), outputBest));
SSres = sum((Dchk(:, end) - outputBest).^2);
SStot = sum((Dchk(:, end) - mean(Dchk(:, end))).^2);
r2 = 1 - SSres/SStot;
nmse = SSres/SStot;
ndei = sqrt(nmse);

% Save numerical results to file
fid = fopen('results/numerical_results.txt', 'w');
fprintf(fid, "----------------\n");
fprintf(fid, "Best TSK Model for radius = %.1f and features = %d\n", bestRadius, bestFeatures);
fprintf(fid, "RMSE = %.4f\n", rmse);
fprintf(fid, "NMSE = %.4f\n", nmse);
fprintf(fid, "NDEI = %.4f\n", ndei);
fprintf(fid, "R^2 = %.4f\n", r2);
fclose(fid);

% Display results in command window
disp("----------------")
fprintf("Best TSK Model for radius = %.1f and features = %d\n", bestRadius, bestFeatures);
fprintf("RMSE = %.4f\n", rmse);
fprintf("NMSE = %.4f\n", nmse);
fprintf("NDEI = %.4f\n", ndei);
fprintf("R^2 = %.4f\n", r2);
