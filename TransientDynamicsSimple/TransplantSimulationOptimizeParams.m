% TransplantSimulationOptimizeParams

% Fit model to data from bare plots in disturbance experiment

clear all

rmpath(genpath(pwd)); % Remove the current working directory and all its subfolders from the path

%% Set path here: %%
path = 'C:\Users\torbe\Nextcloud\Shared\PhD Project\GitCloud\Torben\Chapter TransientDynamics\TransientDynamics\';
%%%%%%%%%%%%%%%%%%%%


%%
delimiter = '';
cd(join([path, 'TransientDynamicsSimple'],delimiter))
%cd('F:\UniCloud\Shared\PhD Project\GitCloud\Torben\Chapter TransientDynamics\TransientDynamics\TransientDynamicsSimple')
addpath(genpath(pwd)); % Add current directory and all subfolders to the path

% IDE Transplant Simulation Model - Optimize parameters
LambdaRgzn = 0.5; % Regularization strength
Species = 4; % Number of species (FG) in the model
randomInitializer = 'n'; % Option to initialize variables randomly within range limits

p = struct();
p.clusters = 4;               % number of clusters to divide community
p.test = 'n';
p.Opt = 'simple';     % Full optimization (with full comp matrix and species specific sigma) % Simple as alternative % Sigma with only Sigma as species specific (sigma not working yet)
p.denConv = 'y';
p.plotInvasion = 'y';

% Parameter templates: Name, initial value, min, max, species-specific flag
paramTemplates = {
    'betaDens',     0.5,     0.01, 10,   false;
    'alphaIntra',   1,       0.1,    2,     false;
    'alphaInter',   1,       0.1,    2,     false;
    'R',            3,       1,      10,    true;  
    'T',            0.5,     0.001,  2,    true;  
    'sigma',        0.05,    0.008,  0.2,   true;
    'gamma',      0.5,       0,      2,     false;
};

% Generate parameter list
paramsTable = createParameterListFromTemplate(paramTemplates, Species);

if randomInitializer == 'y'
    % Replace fixed initial values with random values, within range limits
    paramsTable.InitialValue = arrayfun(@(Min, Max) ...
    Min + rand() * (Max - Min), ...
    paramsTable.Min, paramsTable.Max);
end

% Extract initial values, minimums, and maximums for optimization
initialValues = paramsTable.InitialValue';
paramMin = paramsTable.Min';
paramMax = paramsTable.Max';
VarP_Names = paramsTable.Name;

%% Fixed parameters %%

% Main parameters
p.l = 1;                      % length of the simulated area (from -l/2 to l/2)
p.L = 2^11;                   % Spacial lattice size
p.t_iter = 7;                 % number of iterations (model runtime)
dataToOptimize = 2;
p.Species = p.clusters;       % number of species (FGs)

p.test = 'n';

% Vector with parameters
PVect = paramsTable.InitialValue;

% Optimize cost function
FindMinCost = @(PVect) TransplantSimulationCostFunction(PVect, p, paramsTable, paramTemplates, LambdaRgzn, p.t_iter, dataToOptimize, path);
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 10e+5);  
[FitParam, Delta] = fmincon(FindMinCost, initialValues,[],[],[],[],paramsTable.Min, paramsTable.Max, [], options);     

% After fitting: Test with bare plot data
p.test = 'y';
p.t_iterTest = 7; % Number of years to compare data to
p.simulationName = 'xx'; % Name of simulation, which determines folder name of results
TransplantSimulationCostFunction(FitParam, p, paramsTable, paramTemplates,LambdaRgzn, p.t_iter, dataToOptimize)

% Add fitted values to Parameter table
paramsTable.FitValue = FitParam';


function paramsTable = createParameterListFromTemplate(paramTemplates, Species)
    % Initialize cell array to store parameters
    params = {};
    
    % Loop through the parameter templates
    for i = 1:size(paramTemplates, 1)
        paramName = paramTemplates{i, 1};
        initialValue = paramTemplates{i, 2};
        minVal = paramTemplates{i, 3};
        maxVal = paramTemplates{i, 4};
        isSpeciesSpecific = paramTemplates{i, 5};
        
        if isSpeciesSpecific
            if strcmp(paramName, 'alphaIntra')
                % Special case: create N^2 parameters for 'alphaIntra'
                for sp1 = 1:Species
                    for sp2 = 1:Species
                        % Generate species interaction parameter name (e.g., alphaIntra1_2)
                        newName = sprintf('%s%d_%d', paramName, sp1, sp2);
                        
                        % Base name remains as 'alphaIntra'
                        baseName = paramName;
                        
                        % Append to the params list
                        params = [params; {newName, initialValue, minVal, maxVal, baseName}];
                    end
                end
            else
                % General case: create N parameters for species-specific parameter
                for sp = 1:Species
                    % Generate species-specific parameter name (e.g., R1, R2)
                    newName = sprintf('%s%d', paramName, sp);
                    
                    % Base name is the original parameter name
                    baseName = paramName;
                    
                    % Append to the params list
                    params = [params; {newName, initialValue, minVal, maxVal, baseName}];
                end
            end
        else
            % Add global parameter with base name (same as parameter name)
            baseName = paramName;
            params = [params; {paramName, initialValue, minVal, maxVal, baseName}];
        end
    end
    
    % Convert to table with an additional 'BaseName' column
    paramsTable = cell2table(params, 'VariableNames', {'Name', 'InitialValue', 'Min', 'Max', 'BaseName'});
end
