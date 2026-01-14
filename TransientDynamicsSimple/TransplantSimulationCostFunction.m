function cost = TransplantSimulationCostFunction(PVect, p, paramsTable, paramTemplates, LambdaRgzn, datapoints, dataToOptimize, path)

generateDataStruct = 0; % 1: Option to extract data from PCA table; 0: load data normally
paramsTable.FitValue = PVect'; % Set new values of solver

% Loop through the paramsTable
for i = 1:height(paramsTable)
    % Get the base name and the parameter name
    baseName = paramsTable.BaseName{i};
    paramName = paramsTable.Name{i};
    FitValue = paramsTable.FitValue(i);
    
    % Check if the baseName already exists in the struct p
    if isfield(p, baseName)
        % If it's already in the structure, append the value
        p.(baseName) = [p.(baseName), FitValue];
    else
        % Otherwise, create a new field and assign the value
        p.(baseName) = FitValue;
    end
end


%% Load data

clusterNames = ["Spa", "midM", "ArtFes", "Ely"];
if p.clusters == 3
    clusterNames = ["Spa", "midM", "Ely"];
end
p.clusterNames = clusterNames;

datayears = 7;  % Total number of years in data
dataResCtrl = 200; % Resolution of data along PCA axis (ctrl data)
dataRes = 100;
delimiter = '';
cd(join([path, 'BareData'], delimiter))

%cd('F:\UniCloud\Shared\PhD Project\GitCloud\Torben\Chapter TransientDynamics\TransientDynamics\BareData')

if generateDataStruct == 1 % Daten aus PCA Tabelle extrahieren
    ct = struct(); % Struct where control data is saved
    sm = struct(); % Struct where saltmarsh data is saved
    il = struct(); % Struct where island data is saved
    
    % Should be same size as p.clusterWave
    idx = 0; % Index clusterNr
    ct.x = NaN(datayears,dataResCtrl); % Yr, length(x) % saves x values
    
    for clusterName = clusterNames
        idx = idx +1;
        ct.(clusterName) = NaN(datayears, dataResCtrl); % Yr, length(x) % Saves predictions
        sm.(clusterName) = NaN(datayears, dataRes);
        il.(clusterName) = NaN(datayears, dataRes);
        datact = readtable(sprintf("pred_%s_All_control.csv", clusterName));
        datasm = readtable(sprintf("pred_%s_All_saltmarsh.csv", clusterName));
        if ~strcmp(clusterName, "Ely")
            datail = readtable(sprintf("pred_%s_All_island.csv", clusterName));
        end
        yearidx = 0; % Saves year index
        for year = 2015:2021
            yearidx = yearidx +1;
            datayrct = datact(datact.year == year, :);
            datayrsm = datasm(datasm.year == year, :);
            datayril = datail(datail.year == year, :);
            
            if ~isempty(datayrct)
                ct.(clusterName)(yearidx,:) = datayrct.predicted;
            end
            if ~isempty(datayrsm)
                sm.(clusterName)(yearidx,:) = datayrsm.predicted;
            end
            if ~isempty(datayril)
                il.(clusterName)(yearidx,:) = datayril.predicted;
            end
            
            if idx == 1
                ct.x(yearidx,:) = datayrct.x;
                sm.x(yearidx,:) = datayrsm.x;
                il.x(yearidx,:) = datayril.x;
            end
        end
    end
    il.("Ely")(:) = NaN; % Elymus on island fully to NaN
    save('pcaDataBare200ctrl.mat', 'sm', 'il', 'ct');
else
        load('pcaDataBare200ctrl.mat')
end
smOrig = sm;

% Change working directory back to model
cd(join([path, 'TransientDynamicsSimple'], delimiter))

%cd('C:\Users\torbe\Nextcloud\Shared\PhD Project\GitCloud\Torben\Chapter TransientDynamics\TransientDynamics\TransientDynamicsSimple')
%cd('F:\UniCloud\Shared\PhD Project\GitCloud\Torben\Chapter TransientDynamics\TransientDynamics\TransientDynamicsSimple')


%% Cut data if needed %%

structNames = {'ct', 'sm', 'il'};  % List of struct names as strings

% Cut data if number of datapoints is less than available data
if datapoints < datayears
    for i = 1:length(structNames)
        structname = structNames{i};
        eval([structname '.x = ' structname '.x(1:datapoints, :);']);
        eval([structname '.Spa = ' structname '.Spa(1:datapoints, :);']);
        eval([structname '.midM = ' structname '.midM(1:datapoints, :);']);
        eval([structname '.ArtFes = ' structname '.ArtFes(1:datapoints, :);']);
        eval([structname '.Ely = ' structname '.Ely(1:datapoints, :);']);
    end
end

%% Run model
    p.ct = sm;
    p.control = ct;
    p = TransplantSimulationRunModelCtrlFit(p);

p.densityResults = p.clusterWave; % Save model density results
rawModelResults = p.clusterWave;

    p.clusterWave = p.clusterWave(:,p.ExpModStart:p.ExtHosStart-1,1:p.t_iter); % Only use inside part (without PAD) + dont consider buffertime
    if p.test == 'y'
        p.densityResults = p.densityResults(:,p.ExpModStart:p.ExtHosStart-1,:);
    end


%% select Model results / calculate frequency

% Estimate frequency based on density
p.clusterWave = 1-exp(-p.betaDens'.*p.clusterWave);
if p.test == 'y'
    p.densityResults = p.clusterWave;
end

%% Remove model data where actual data is NaN
m = struct(); % Struct for model results
mTest = struct();

structIdx = 0; % Index for structures
for i = 1:length(structNames)  % loop over ct, sm, il
    structIdx = structIdx + 1;  
    structname = structNames{i};
    clusterIdx = 0;
    for clusterName = clusterNames % loop over all clusters
        clusterIdx = clusterIdx + 1;
        
        rowsWithNaN = any(isnan(eval([structname '.(clusterName)'])), 2); % rows in data that contain NaN (for specific cluster)
        eval([structname '.(clusterName) =' structname '.(clusterName)(~rowsWithNaN, :);']) % only select data without NaN's
        
        if structIdx == dataToOptimize % if cluster is considered that should be used for cost calculation
            m.(clusterName) = squeeze(p.clusterWave(clusterIdx,:,~rowsWithNaN))'; % save model results of cluster without rows where data has NaN
        end
    end
end

if p.test == 'y'
    structname = 'smOrig';
    clusterIdx = 0;
    for clusterName = clusterNames % loop over all clusters
        clusterIdx = clusterIdx + 1;
        
        rowsWithNaN = any(isnan(eval([structname '.(clusterName)'])), 2); % rows in data that contain NaN (for specific cluster)
        rowsWithNaN = rowsWithNaN(1:datapoints); % Only look at range of datapoints that shall be considered
        eval([structname '.(clusterName) =' structname '.(clusterName)(~rowsWithNaN, :);']) % only select data without NaN's
        
        mTest.(clusterName) = squeeze(p.densityResults(clusterIdx,:,~rowsWithNaN))'; % save model results of cluster without rows where data has NaN
        if p.test == 'y'
            mTest.(clusterName) = squeeze(p.densityResults(clusterIdx,:,~rowsWithNaN))'; % save model results of cluster without rows where data has NaN
        end
    end
end

%% Rescale spatial size of model to data
modelResolution = p.L*(2/8);
dataResolution = size(sm.Spa,2);

xNew = linspace(1, modelResolution, dataResolution);
xOriginal = linspace(1, modelResolution, modelResolution);

% Interpolate each row
for clusterName = clusterNames % Loop for each cluster
    m.New = zeros(size(m.(clusterName),1), dataResolution); % Initialize matrix
    if p.test == 'y'
        mTest.New = zeros(size(mTest.(clusterName),1), dataResolution); % Initialize matrix
    end
    
    for i = 1:size(m.(clusterName),1) % Loop for each year
        m.New(i,:) = interp1(xOriginal, m.(clusterName)(i,:), xNew);
    end
    if p.test == 'y'
        for i = 1:size(smOrig.Spa,1) % Loop for each year
            mTest.New(i,:) = interp1(xOriginal, mTest.(clusterName)(i,:), xNew);
        end
    end
    m.(clusterName) = m.New;
    if p.test == 'y'
        mTest.(clusterName) = mTest.New;
    end
end

%% Invert data --> Benign left, Stress right
dataStruct = eval([structNames{dataToOptimize}]); % Set significant datastruct depending on which type should be optimized
for clusterName = clusterNames % Loop for each cluster
    dataStruct.(clusterName) =  fliplr(dataStruct.(clusterName));
end
testDataStruct = smOrig;

%% Calculate cost
costSum = struct();
for clusterName = clusterNames % Loop over all clusters
    costSum.(clusterName) = sum((m.(clusterName)-dataStruct.(clusterName)).^2,'all');
    costSum.(sprintf('Year%s', clusterName)) = sum((m.(clusterName)-dataStruct.(clusterName)).^2,2);
end

% With Regularization
RegCost = 0; % Set regularization cost initially to 0
for i = 1:length(paramsTable.Name)
    for s = 1:size(p.R)
        p.(sprintf("R%i", s)) = p.R(s);
        p.(sprintf("T%i", s)) = p.T(s);
        p.(sprintf("sigma%i", s)) = p.sig(s);
        for j = 1:size(p.R)
            p.(sprintf("alpha%i%i", s,j)) = p.comp(s,j);
        end
    end
    meanParam = mean([paramsTable.Min(i), paramsTable.Max(i)]);
    deltaParam = abs(paramsTable.Max(i)-paramsTable.Min(i));
    RegCost = RegCost + LambdaRgzn * (abs(paramsTable.FitValue(i) - meanParam) ./ deltaParam);
end
costSum.RegCost = RegCost;

% Square cost of each cluster
costSum.YearSum = zeros(p.t_iter,1);
for clusterName = clusterNames
    costSum.YearSum = costSum.YearSum + costSum.(sprintf('Year%s', clusterName));
end
cost = sum(costSum.YearSum); % Combined cost of difference between model and data
cost = cost + RegCost; % Add Regularization cost

%% Generate comparison plot

if p.test == 'y'
    p.t_iterTest = size(smOrig.Spa,1);
    p.t_iter = p.t_iterTest;

    % Frequencies over time
    f_t = figure();
    load('FrequencyData.mat')
    hold on
    samplePoints = [0.1,0.5,0.9];
    sampleIdx = round(size(mTest.Spa,2).*samplePoints);
    for clusterName = clusterNames
        for zone = 1:length(samplePoints)
            idx = sampleIdx(zone);
            subplot(2,length(samplePoints),zone)
            plot(1:datapoints,mTest.(clusterName)(:,idx), 'LineWidth', 2)
            hold on
            ylim([0,1])
            xticks([0,1,2,3,4,5,6,7])
            %xlabel('t', 'FontSize', 14)
            if zone == 1
                title('Upper', 'FontSize', 12)
                ylabel('Model f_{m,i}', 'FontSize',14)
            elseif zone == 2
                title('Lower', 'FontSize', 12)
            elseif zone == 3
                title('Pio', 'FontSize', 12)
            end
            subplot(2, length(samplePoints), 2*length(samplePoints)+1-zone)
            tab = BareAllElevations{zone}; % Select table of zone
            plot(1:datapoints, tab.(clusterName)(1:datapoints), 'LineWidth', 2)
            xlabel('t', 'FontSize', 14)
            if zone == 3
                ylabel('Data f_{d,i}', 'FontSize', 14)
            end
            xticks([0,1,2,3,4,5,6,7])
            hold on
        end
    end
    leg = legend(clusterNames, 'FontSize', 10);
    leg.ItemTokenSize = [15, 10];
    hold off
    AddLetters2Plots(gcf, {'a)', 'b)', 'c)', 'd)','e)', 'f)'}, 'HShift', -0.065, 'VShift', -0.065)
    
    
    % Comparison Model vs Data - Distributions across gradient
    % Generate a colormap for the iterations
    
    f_comp = figure();
    %colormapYears = lines(p.t_iter); % Colormap for both model and empirical data
    colormapYears = [
    0.231, 0.298, 0.753;  % deep blue
    0.344, 0.533, 0.949;  % medium blue
    0.196, 0.698, 0.494;  % teal
    0.596, 0.875, 0.541;  % light green
    0.996, 0.835, 0.369;  % yellow
    0.992, 0.482, 0.223;  % orange
    0.843, 0.153, 0.122;  % red
];

    
    sidx = length(clusterNames) + 1; % Subplot index
    for clusterName = clusterNames
        sidx = sidx - 1;
        subplot(2, length(clusterNames)/2, sidx)
        hold on
        % Pre-allocate legend entries for years
        legendEntries = cell(1, p.t_iter);
        dummyLines = gobjects(1, p.t_iter); % Dummy lines for consistent legend styles
        for t = 1:p.t_iter
            % Plot model data with dashed lines
            plot(linspace(0, 1), mTest.(clusterName)(t, :), '--', 'Color', colormapYears(t, :), 'LineWidth', 1.5)
            % Plot empirical data with solid lines
            plot(linspace(0, 1), fliplr(testDataStruct.(clusterName)(t, :)), '-', 'Color', colormapYears(t, :), 'LineWidth', 1.5)
            % Add dummy line for consistent year legend (always solid style)
            dummyLines(t) = plot(NaN, NaN, '-', 'Color', colormapYears(t, :), 'LineWidth', 1.5);
            % Add legend entry for the current year
            legendEntries{t} = sprintf('Year %d', t);
        end
        title(sprintf('FG %s', clusterName), 'FontSize', 12)
        ylabel('Freq. f_{i,t}(x)', 'FontSize', 14)
        xlabel('x', 'FontSize', 14)
        if sidx == 1
            % add legend
            legend(dummyLines, legendEntries, 'Location', 'best', 'FontSize', 10)
        end
        ylim([0,1])
        if strcmp(clusterName, 'ArtFes')
            ylim([0, 0.4])
        end
        hold off        
    end
    figureHandle = gcf; % Get the current figure
    figureHandle.Position = [70,24,826,580];
    
    % Add a global legend for line styles
    axesHandle = axes('Position', [0, 0, 1, 1], 'Visible', 'off'); % Invisible axes for a global legend
    hold on
    plot(NaN, NaN, '--k', 'LineWidth', 1.5); % Example line for model
    plot(NaN, NaN, '-k', 'LineWidth', 1.5);  % Example line for empirical data
    legend(axesHandle, {'Model (dashed)', 'Empirical (solid)'}, 'Location', 'southoutside', 'Orientation', 'vertical', 'FontSize', 10)
    hold off
    AddLetters2Plots(figureHandle, {'a)', 'b)', 'c)', 'd)'}, 'HShift', -0.06, 'VShift', -0.05, 'FontSize', 14)

    %sgtitle('Model vs. Test Data')
    drawnow
    
    % Save Figures and ParamTable
    saveResults(f_t, f_comp, paramsTable, p.simulationName)
       
    % Generate results and simulationName folder (if not already existing)
    % Ensure the Parent 'Results' Directory Exists
    parentFolder = fullfile(pwd, 'Results'); % Path to the Results folder
    if ~exist(parentFolder, 'dir')
        mkdir(parentFolder);
    end
    
    % Ensure the 'simulationName' Subfolder Exists
    simulationFolder = fullfile(parentFolder, p.simulationName);
    if ~exist(simulationFolder, 'dir')
        mkdir(simulationFolder);
    end
    
    %Save workspace
    workspaceFile = fullfile(simulationFolder, sprintf('%s.mat', p.simulationName));
    save(workspaceFile);
    
end

end

