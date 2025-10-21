close all
clear
clc

% Load Data and Start Preprocessing
data = importdata("epileptic_seizure_data.csv");
data = data.data;
[m, n] = size(data);

preproc = 1;
[Dtrn, Dval, Dchk] = stratified_split(data, preproc);

% Define kFold, Radius and Features
KFold = 5;
feature = [5 8 11 15];
radius = [0.3 0.5 0.7 0.9];

n1 = length(feature);
n2 = length(radius);

[indexList, weights] = relieff(Dtrn(:, 1:end - 1), Dtrn(:, end), 10);

for i = 1: 1: n1
    for j = 1: 1: n2
        cvpart = cvpartition(length(Dtrn(:, end)), 'KFold', KFold);
        for k = 1: 1: KFold
            fprintf("radius = %d, feature = %d", radius(j), feature(i));
            % Choose Features
            chosen = [Dtrn(:, indexList(1:feature(i))) Dtrn(:, end)];

            index1 = test(cvpart, k); 
            index2 = training(cvpart, k);

            validationIndex = find(index1 == 1);
            trainIndex = find(index2 == 1);

            DtrnNew = chosen(trainIndex, :);
            DvalNew = chosen(validationIndex, :);

            [rows, columns] = size(DtrnNew);

            % Find Cluster Centers using Subtractive Clustering
            [c1, sig1] = subclust(DtrnNew(DtrnNew(:, end) == 1, :), radius(j));
            [c2, sig2] = subclust(DtrnNew(DtrnNew(:, end) == 2, :), radius(j));
            [c3, sig3] = subclust(DtrnNew(DtrnNew(:, end) == 3, :), radius(j));
            [c4, sig4] = subclust(DtrnNew(DtrnNew(:, end) == 4, :), radius(j));
            [c5, sig5] = subclust(DtrnNew(DtrnNew(:, end) == 5, :), radius(j));
            num_rules = size(c1, 1) + size(c2, 1) + size(c3, 1) + size(c4, 1) + size(c5, 1);

            % Create Fuzzy Inference System (FIS)
            tsk = sugfis;

            % Add Input Variables
            for ii = 1: 1: columns - 1
                nameGiven = "In" + int2str(ii);
                tsk = addInput(tsk, [0 1], 'Name', nameGiven);
                nameGivenList(ii) = nameGiven;
            end

            % Add Output Variables
            tsk = addOutput(tsk, [1 5], 'Name', 'Output');

            % Add Input Membership Functions

            for ii = 1: 1: columns - 1
                for jj = 1: 1: size(c1, 1)
                    tsk = addMF(tsk, nameGivenList(ii), 'gaussmf', [sig1(ii) c1(jj, ii)]);
                end

                for jj = 1: 1: size(c2, 1)
                    tsk = addMF(tsk, nameGivenList(ii), 'gaussmf', [sig2(ii) c2(jj, ii)]);
                end

                for jj = 1: 1: size(c3, 1)
                    tsk = addMF(tsk, nameGivenList(ii), 'gaussmf', [sig3(ii) c3(jj, ii)]);
                end

                for jj = 1: 1: size(c4, 1)
                    tsk = addMF(tsk, nameGivenList(ii), 'gaussmf', [sig4(ii) c4(jj, ii)]);
                end

                for jj = 1: 1: size(c5, 1)
                    tsk = addMF(tsk, nameGivenList(ii), 'gaussmf', [sig5(ii) c5(jj, ii)]);
                end
            end

            % Add Output Membership Functions
            parameters = [ones(1, size(c1, 1)) 2 * ones(1, size(c2, 1)) 3 * ones(1, size(c3, 1)) 4 * ones(1, size(c4, 1)) 5 * ones(1, size(c5, 1))];

            for ii = 1: 1: num_rules
                tsk = addMF(tsk, 'Output', 'constant', parameters(ii));
            end

            % Adding the Rules
            ruleList = zeros(num_rules, size(DtrnNew, 2));

            for ii = 1: 1: size(ruleList, 1)
                ruleList(ii, :) = ii;
            end
            ruleList = [ruleList ones(num_rules, 2)];
            tsk = addrule(tsk, ruleList);

            % Train Fuzzy Inference System (FIS)
            options2 = anfisOptions;
            options2.InitialFIS = tsk;
            options2.EpochNumber = 100;
            options2.ValidationData = DvalNew;

            [fis, trainError, stepSize, valFIS, valError] = anfis(DtrnNew, options2);

            rules(k) = size(showrule(valFIS), 1);
            errors(k) = abs(mean(valError));
        end
        errors_(i, j) = mean(errors);
        rules_(j) = round(mean(rules));
    end    
end


% Plots Mean Absolute Error - Rules and Mean Absolute Error - Features
hold on
grid on
title('Mean Error - Number Of Features', 'Interpreter', 'latex');
xlabel('Number Of Features', 'Interpreter', 'latex')
ylabel('Mean Absolute Error', 'Interpreter', 'latex')
plot(feature, errors_(:, 1), 'o', 'LineWidth', 1, 'Color', 'red');
plot(feature, errors_(:, 2), 'o', 'LineWidth', 1, 'Color', 'blue');
plot(feature, errors_(:, 3), 'o', 'LineWidth', 1, 'Color', 'green');
plot(feature, errors_(:, 4), 'o', 'LineWidth', 1, 'Color', 'black');
labels = num2str([rules_(1), rules_(2), rules_(3), rules_(4)].', '%d Rules');
legend(labels);
figure

hold on
grid on
title('Mean Error - Number Of Rules', 'Interpreter', 'latex');
xlabel('Number Of Rules', 'Interpreter', 'latex');
ylabel('Mean Error', 'Interpreter', 'latex');
plot(rules_, errors_(1, :), 'o', 'LineWidth', 1, 'Color', 'red');
plot(rules_, errors_(2, :), 'o', 'LineWidth', 1, 'Color', 'blue');
plot(rules_, errors_(3, :), 'o', 'LineWidth', 1, 'Color', 'green');
plot(rules_, errors_(4, :), 'o', 'LineWidth', 1, 'Color', 'black');
labels = num2str([feature(1), feature(2), feature(3), feature(4)].', '%d Features');
legend(labels);
figure

% Get Best Model
[bestIndex1, bestIndex2] = find(errors_ == min(min(errors_)));

bestFeatures = feature(bestIndex1);
bestRadius = radius(bestIndex2);

DtrnBest = [Dtrn(:, indexList(1:bestFeatures)) Dtrn(:, end)];
DvalBest = [Dval(:, indexList(1:bestFeatures)) Dval(:, end)];
DchkBest = [Dchk(:, indexList(1:bestFeatures)) Dchk(:, end)];

[c1best, sig1best] = subclust(DtrnBest(DtrnBest(:, end) == 1, :), bestRadius);
[c2best, sig2best] = subclust(DtrnBest(DtrnBest(:, end) == 2, :), bestRadius);
[c3best, sig3best] = subclust(DtrnBest(DtrnBest(:, end) == 3, :), bestRadius);
[c4best, sig4best] = subclust(DtrnBest(DtrnBest(:, end) == 4, :), bestRadius);
[c5best, sig5best] = subclust(DtrnBest(DtrnBest(:, end) == 5, :), bestRadius);
best_num_rules = size(c1best, 1) + size(c2best, 1) + size(c3best, 1) + size(c4best, 1) + size(c5best, 1);

% Create Fuzzy Inference System(FIS)
bestTsk = sugfis;

% Add Input-Output Variables
[rowsBest, columnsBest] = size(DtrnBest);

for ii = 1: 1: columnsBest - 1
    nameGivenBest = "In" + int2str(ii);
    bestTsk = addInput(bestTsk, [0 1], 'Name', nameGivenBest);
    nameGivenListBest(ii) = nameGivenBest;
end
bestTsk = addOutput(bestTsk, [1 5], 'Name', 'OutputBest');

% Add Input Membership Functions
for ii = 1: 1: columnsBest - 1
    for jj = 1: 1: size(c1best,1)
        bestTsk = addMF(bestTsk, nameGivenListBest(ii), 'gaussmf', [sig1best(ii) c1best(jj,ii)]);
    end
    for jj = 1: 1: size(c2best, 1)
        bestTsk = addMF(bestTsk, nameGivenListBest(ii), 'gaussmf', [sig2best(ii) c2best(jj,ii)]);
    end
    for jj = 1: 1: size(c3best, 1)
        bestTsk = addMF(bestTsk, nameGivenListBest(ii), 'gaussmf', [sig3best(ii) c3best(jj,ii)]);
    end
    for jj = 1: 1: size(c4best, 1)
        bestTsk = addMF(bestTsk, nameGivenListBest(ii), 'gaussmf', [sig4best(ii) c4best(jj,ii)]);
    end
    for jj = 1: 1: size(c5best, 1)
        bestTsk = addMF(bestTsk, nameGivenListBest(ii), 'gaussmf', [sig5best(ii) c5best(jj,ii)]);
    end
end

% Add Output Membership Functions
parametersBest = [ones(1, size(c1best, 1)) 2 * ones(1, size(c2best, 1)) 3 * ones(1, size(c3best, 1)) 4 * ones(1, size(c4best, 1)) 5 * ones(1, size(c5best, 1))];
for ii = 1: 1: best_num_rules
    bestTsk = addMF(bestTsk,'OutputBest', 'constant', parametersBest(ii));
end

% Adding the Rules
best_ruleList = zeros(best_num_rules, columnsBest);
for ii = 1: 1: size(best_ruleList, 1)
    best_ruleList(ii, :) = ii;
end
best_ruleList = [best_ruleList ones(best_num_rules, 2)];
bestTsk = addrule(bestTsk, best_ruleList);

% Train Fuzzy Inference System
optionsBest2 = anfisOptions;
optionsBest2.InitialFIS = bestTsk;
optionsBest2.EpochNumber = 100;
optionsBest2.ValidationData = DvalBest;

[fis, trainErrorB, stepSize, valFIS, valErrorB] = anfis(DtrnBest, optionsBest2);

% Evaluate Fuzzy Inference System
outputBest = round(evalfis(DchkBest(:, 1:end - 1), valFIS));

% Calculate ErrorMatrix, OverallAccuracy, ProducerAccuracy, UserAccuracy, Khat
classes = [1; 2; 3; 4; 5];
ErrorMatrix = zeros(5, 5);

for ii = 1: 1: length(DchkBest)
    a = find(classes == outputBest(ii));
    b = find(classes == data(ii, end));
    ErrorMatrix(a, b) = ErrorMatrix(a, b) + 1;
end

N = length(DchkBest);

OA = 1/N * sum(diag(ErrorMatrix));

sum1 = sum(ErrorMatrix, 2);
sum2 = sum(ErrorMatrix, 1);

for i = 1: 1: 5
    PA(i) = ErrorMatrix(i, i) / sum2(i);
    UA(i) = ErrorMatrix(i, i) / sum1(i);
end

kHat = (N * trace(ErrorMatrix) - PA * UA')/(N^2 - PA * UA');


totalNumberOfRules = size(valFIS.Rules,2);

figure
hold on
grid on
title('Mean Error - Number Of Features', 'Interpreter', 'latex');
xlabel('Number Of Features', 'Interpreter', 'latex')
ylabel('Mean Absolute Error', 'Interpreter', 'latex')
plot(feature, errors_(:, 1), 'o-', 'LineWidth', 1, 'Color', 'red');
plot(feature, errors_(:, 2), 'o-', 'LineWidth', 1, 'Color', 'blue');
plot(feature, errors_(:, 3), 'o-', 'LineWidth', 1, 'Color', 'green');
plot(feature, errors_(:, 4), 'o-', 'LineWidth', 1, 'Color', 'black');
legend('radius=0.3', 'radius=0.5', 'radius=0.7', 'radius=0.9');

figure
hold on
grid on
title('Mean Error - Number Of Rules', 'Interpreter', 'latex');
xlabel('Number Of Rules', 'Interpreter', 'latex');
ylabel('Mean Error', 'Interpreter', 'latex');
plot(rules_, errors_(1, :), 'o-', 'LineWidth', 1, 'Color', 'red');
plot(rules_, errors_(2, :), 'o-', 'LineWidth', 1, 'Color', 'blue');
plot(rules_, errors_(3, :), 'o-', 'LineWidth', 1, 'Color', 'green');
plot(rules_, errors_(4, :), 'o-', 'LineWidth', 1, 'Color', 'black');
legend('5 features', '8 features', '11 features', '15 features');

for i = 1:2
    figure
    subplot(2,1,1)
    plotmf(bestTsk, 'input', i);
    title("Before Training", 'interpreter','latex');
    subplot(2,1,2)
    plotmf(valFIS, 'input', i);
    title("After Training", 'interpreter','latex');
end

figure
plot([trainErrorB valErrorB]);
legend('Training Error', 'Validation Error');
title("Learning Curve", 'interpreter','latex');

figure
plot(1:length(outputBest), outputBest, 'r*');
title("Predicted Values", 'interpreter','latex');

figure
plot(1:length(outputBest), DchkBest(:, end), 'bo');
title("Real Values", 'interpreter','latex');

figure
plot(DchkBest(:, end) - outputBest);
title("Prediction Error", 'interpreter','latex');

% Display ErrorMatrix, OverallAccuracy, ProducerAccuracy, UserAccuracy, Khat
disp("--------------------------------------");
fprintf("Best TSK Model for total number of features = %d and cluser radius = %d\n", bestFeatures, bestRadius);
fprintf("Total Number Of Rules is = %d\n", totalNumberOfRules);
fprintf("OA = %d\n", OA);
fprintf("PA = %d\n", PA);
fprintf("UA = %d\n", UA);
fprintf("KHat = %d\n", kHat);
disp("Error Matrix");
disp(ErrorMatrix);
