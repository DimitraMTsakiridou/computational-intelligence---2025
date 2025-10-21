close all; clear variables; clc;

% Directory setup for saving outputs
plotDirectory = 'model_plots';
if ~isfolder(plotDirectory)
    mkdir(plotDirectory);
end

% Import and prepare dataset
dataset = importdata("airfoil_self_noise.dat");
[rows, cols] = size(dataset);

preprocessingMode = 1;
[TrainingData, ValidationData, TestingData] = split_scale(dataset, preprocessingMode);

% Model configuration parameters
membershipFunctionTypes = ["constant", "constant", "linear", "linear"];
numMFs = [2, 3, 2, 3];
totalModels = 4;

% Initialize performance metrics storage
performanceMetrics = struct();
metricNames = {'rmse', 'nmse', 'ndei', 'r2'};
for m = 1:length(metricNames)
    performanceMetrics.(metricNames{m}) = zeros(1, totalModels);
end

% Main modeling loop
for modelIdx = 1:totalModels
    % FIS generation configuration
    fisOptions = genfisOptions('GridPartition');
    fisOptions.NumMembershipFunctions = numMFs(modelIdx);
    fisOptions.InputMembershipFunctionType = 'gbellmf';
    fisOptions.OutputMembershipFunctionType = membershipFunctionTypes(modelIdx);
    
    initialFIS = genfis(TrainingData(:, 1:end-1), TrainingData(:, end), fisOptions);
    
    % Training configuration
    trainingOpts = anfisOptions;
    trainingOpts.InitialFIS = initialFIS;
    trainingOpts.ValidationData = ValidationData;
    trainingOpts.EpochNumber = 100;  
    trainingOpts.OptimizationMethod = 1;
    
    [trainedFIS, trainingError, ~, valFIS, validationError] = anfis(TrainingData, trainingOpts);
    
    % Model evaluation
    predictions = evalfis(valFIS, TestingData(:, 1:end-1));
    actualValues = TestingData(:, end);
    
    % Calculate performance metrics
    performanceMetrics.rmse(modelIdx) = sqrt(mean((actualValues - predictions).^2));
    ss_residual = sum((actualValues - predictions).^2);
    ss_total = sum((actualValues - mean(actualValues)).^2);
    performanceMetrics.r2(modelIdx) = 1 - (ss_residual/ss_total);
    performanceMetrics.nmse(modelIdx) = ss_residual/ss_total;
    performanceMetrics.ndei(modelIdx) = sqrt(performanceMetrics.nmse(modelIdx));
    
    % Visualization of membership functions
    for featureNum = 1:5
        fig = figure('Visible', 'off');
        subplot(2, 1, 1);
        plotmf(initialFIS, 'input', featureNum);
        title(sprintf('Initial MFs for Feature %d', featureNum), 'Interpreter', 'latex');
        
        subplot(2, 1, 2);
        plotmf(trainedFIS, 'input', featureNum);
        title(sprintf('Trained MFs for Feature %d', featureNum), 'Interpreter', 'latex');
        
        saveas(fig, fullfile(plotDirectory, sprintf('model%d_feature%d_membership.png', modelIdx, featureNum)));
        close(fig);
    end
    
    % Training progress visualization
    fig = figure('Visible', 'off');
    plot([trainingError validationError], 'LineWidth', 1.5);
    legend({'Training Error', 'Validation Error'}, 'Location', 'best');
    title('Model Training Progress', 'Interpreter', 'latex');
    xlabel('Epochs');
    ylabel('Error');
    saveas(fig, fullfile(plotDirectory, sprintf('model%d_training_progress.png', modelIdx)));
    close(fig);
    
    % Prediction error analysis
    fig = figure('Visible', 'off');
    errorValues = actualValues - predictions;
    plot(errorValues, 'LineWidth', 1.5);
    title('Model Prediction Errors', 'Interpreter', 'latex');
    xlabel('Sample Index');
    ylabel('Error Value');
    saveas(fig, fullfile(plotDirectory, sprintf('model%d_prediction_errors.png', modelIdx)));
    close(fig);
end

% Display performance results
for k = 1:totalModels
    fprintf('\n--------------------------------\n');
    fprintf('Performance of Model %d:\n', k);
    fprintf('Root Mean Squared Error: %.4f\n', performanceMetrics.rmse(k));
    fprintf('Normalized Mean Squared Error: %.4f\n', performanceMetrics.nmse(k));
    fprintf('Non-Dimensional Error Index: %.4f\n', performanceMetrics.ndei(k));
    fprintf('R-squared Coefficient: %.4f\n', performanceMetrics.r2(k));
end

% Save results to file
resultsFile = fopen(fullfile(plotDirectory, 'performance_results.txt'), 'w');
for k = 1:totalModels
    fprintf(resultsFile, '\n--------------------------------\n');
    fprintf(resultsFile, 'Performance of Model %d:\n', k);
    fprintf(resultsFile, 'Root Mean Squared Error: %.4f\n', performanceMetrics.rmse(k));
    fprintf(resultsFile, 'Normalized Mean Squared Error: %.4f\n', performanceMetrics.nmse(k));
    fprintf(resultsFile, 'Non-Dimensional Error Index: %.4f\n', performanceMetrics.ndei(k));
    fprintf(resultsFile, 'R-squared Coefficient: %.4f\n', performanceMetrics.r2(k));
end
fclose(resultsFile);