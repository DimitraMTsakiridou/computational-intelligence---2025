close all
clear
clc
format LongG

% Load Data and Start Preprocessing
data = importdata('haberman.data');
[m, n] = size(data);

preproc = 1;
[Dtrn, Dval, Dchk] = stratified_split(data, preproc);

% Define Radius
radius = [0.3 0.8];

% Class Dependent
for i = 1: 1: 2

    [rows, columns] = size(Dtrn);

    % Find Cluster Centers using Subtractive Clustering
    [c1, sig1] = subclust(Dtrn(Dtrn(:, end) == 1, :), radius(i));
    [c2, sig2] = subclust(Dtrn(Dtrn(:, end) == 2, :), radius(i));
    numRules = size(c1, 1) + size(c2, 1);

    % Create Sugeno Fuzzy Inference System (FIS)
    tsk = sugfis;

    % Add Input Variables
    for j = 1: 1: columns - 1
        nameGiven = "In" + int2str(j);
        tsk = addInput(tsk, [0 1], 'Name', nameGiven);
        nameGivenList(j) = nameGiven;
    end

    % Add Output Variables
    tsk = addOutput(tsk, [1 2], 'Name', 'Output');

    % Add Input Membership Functions 
    for j = 1: 1: columns - 1
        for k = 1: 1: size(c1, 1)    
            tsk = addMF(tsk, nameGivenList(j), 'gaussmf', [sig1(j) c1(k, j)]);
        end
        for l = 1: 1: size(c2, 1)
            tsk = addMF(tsk, nameGivenList(j), 'gaussmf', [sig2(j) c2(l, j)]);
        end
    end
    
    % Add Output Membership Functions
    parameters = [ones(1, size(c1, 1)) 2 * ones(1, size(c2, 1))];

    for j = 1: 1: numRules
        tsk = addMF(tsk, 'Output', 'constant', parameters(j));
    end

    % Adding the Rules
    ruleList = zeros(numRules, columns);
    for j = 1: 1: size(ruleList, 1)
        ruleList(j, :) = j;
    end
    ruleList = [ruleList ones(numRules, 2)];
    
    tsk = addrule(tsk, ruleList);
  
    % Train Fuzzy Inference System (FIS)
    options2 = anfisOptions;
    options2.InitialFIS = tsk;
    options2.EpochNumber = 100;
    options2.ValidationData = Dval;

    [fis, trainError, stepSize, valFIS, valError] = anfis(Dtrn, options2);

    % Evaluate Fuzzy Inference System (FIS)
    output = round(evalfis(valFIS, Dchk(:, 1:end - 1)));
    output = min(max(1, output), 2);

    % Calculate ErrorMatrix, OverallAccuracy, ProducerAccuracy, UserAccuracy, Khat

    classes = [1; 2];

    ErrorMatrix = zeros(2, 2);

    for ii = 1: 1: length(Dchk)
        a = find(classes == output(ii));
        b = find(classes == data(ii, end));
        ErrorMatrix(a, b) = ErrorMatrix(a, b) + 1;
    end
    
    N = length(Dchk);
    OA = 1/N * trace(ErrorMatrix);

    sum1 = sum(ErrorMatrix, 2);
    sum2 = sum(ErrorMatrix, 1);

    for j = 1: 1: 2
        UA(j) = ErrorMatrix(j, j) / sum1(j);
        PA(j) = ErrorMatrix(j, j) / sum2(j);
    end

    kHat = (N * trace(ErrorMatrix) - PA * UA')/(N^2 - PA * UA');

    OA_array(i) = OA;
    PA_array(i, :, :) = PA;
    UA_array(i, :, :) = UA;
    kHat_array(i) = kHat;
    ErrorMatrixArray(i, :, :) = ErrorMatrix;

    totalNumberOfRules(i) = size(valFIS.Rules, 2);

    % Plot Membership Function after Training
    plotmf(fis, 'input', 1);
    title(sprintf("feature 1 after training - Class Dependent with Radius %d", radius(i)), 'interpreter','latex');
    figure
    plotmf(fis, 'input', 2);
    title(sprintf("feature 2 after training - Class Dependent with Radius %d", radius(i)), 'interpreter','latex');    
    figure
    plotmf(fis, 'input', 3);
    title(sprintf("feature 3 after training - Class Dependent with Radius %d", radius(i)), 'interpreter','latex');
    figure
 
    % Plot Learning Curve
    plot([trainError valError]);
    legend('Training Error', 'Validation Error');
    title("Learning Curve - Class Dependent", 'interpreter','latex');
    figure;

end

% Class Independent
for i = 1: 1: 2

    % Create Sugeno Fuzzy Inference System (FIS)
    options = genfisOptions('SubtractiveClustering');
    options.ClusterInfluenceRange = radius(i);
    tsk = genfis(Dtrn(:, 1:end - 1), Dtrn(:, end), options);

    % Train Fuzzy Inference System (FIS)
    options2 = anfisOptions;
    options2.InitialFIS = tsk;
    options2.EpochNumber = 100;
    options2.ValidationData = Dval;

    [fis, trainError, stepSize, valFIS, valError] = anfis(Dtrn, options2);

    % Evaluate Fuzzy Inference System (FIS)
    output = round(evalfis(valFIS, Dchk(:, 1:end - 1)));
    output = min(max(1, output), 2);

    % Calculate ErrorMatrix, OverallAccuracy, ProducerAccuracy, UserAccuracy, Khat
    classes = [1; 2];

    ErrorMatrix2 = zeros(2, 2);

    for ii = 1: 1: length(Dchk)
        a = find(classes == output(ii));
        b = find(classes == data(ii, end));
        ErrorMatrix2(a, b) = ErrorMatrix2(a, b) + 1;
    end
    
    N2 = length(Dchk);
    OA2 = 1/N2 * sum(diag(ErrorMatrix2));

    sum1_2 = sum(ErrorMatrix2, 2);
    sum2_2 = sum(ErrorMatrix2, 1);

    for j = 1: 1: 2
        UA2(j) = ErrorMatrix2(j, j) / sum1_2(j);
        PA2(j) = ErrorMatrix2(j, j) / sum2_2(j);
    end

    kHat2 = (N2 * trace(ErrorMatrix2) - PA2 * UA2')/(N2^2 - PA2 * UA2');

    OA2_array(i) = OA2;
    PA2_array(i, :, :) = PA2;
    UA2_array(i, :, :) = UA2;
    kHat2_array(i) = kHat2;
    ErrorMatrixArray2(i, :, :) = ErrorMatrix2;

    totalNumberOfRules2(i) = size(valFIS.Rules, 2);

    % Plot Membership Function after Training
    plotmf(fis, 'input', 1);
    title(sprintf("feature 1 after training - Class Independent with Radius %d", radius(i)), 'interpreter','latex');
    figure
    plotmf(fis, 'input', 2);
    title(sprintf("feature 2 after training - Class Independent with Radius %d", radius(i)), 'interpreter','latex');    
    figure
    plotmf(fis, 'input', 3);
    title(sprintf("feature 3 after training - Class Independent with Radius %d", radius(i)), 'interpreter','latex');
    figure

    % Plot Learning Curve
    plot([trainError valError]);
    legend('Training Error', 'Validation Error');
    title("Learning Curve - Class Independent", 'interpreter','latex');
    
    if i ~= 2
        figure
    end
end

% Display ErrorMatrix, OverallAccuracy, ProducerAccuracy, UserAccuracy, Khat
for k = 1: 1: 2
    disp("--------------------------------------");
    fprintf("Class Dependent with radius = %d\n", radius(k));
    fprintf("OA = %d\n", OA_array(k));
    fprintf("PA = %d\n", PA_array(k, :, :));
    fprintf("UA = %d\n", UA_array(k, :, :));
    fprintf("KHat = %d\n", kHat_array(k));
    fprintf("Error Matrix is \n");
    disp(ErrorMatrixArray(:, :, k));
    fprintf("Total Number Of Rules = %d\n", totalNumberOfRules(k));
end
for k = 1: 1: 2
    disp("--------------------------------------");
    fprintf("Class Independent with radius = %d\n", radius(k));
    fprintf("OA = %d\n", OA2_array(k));
    fprintf("PA = %d\n", PA2_array(k, :, :));
    fprintf("UA = %d\n", UA2_array(k, :, :));
    fprintf("KHat = %d\n", kHat2_array(k));
    fprintf("Error Matrix is \n");
    disp(ErrorMatrixArray2(:, :, k));
    fprintf("Total Number Of Rules = %d\n", totalNumberOfRules2(k));
end