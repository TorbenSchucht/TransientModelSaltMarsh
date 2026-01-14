function p = TransplantSimulationRunModelCtrlFit(p)

% IDE Transplant Simulation Model - Runs model with parameters, returns N
% Simulates invasion into open space with species from a Metacommunity

%%%%%%%%%%%%%%%%%%%%
%% DEFINE METHODS %%
%%%%%%%%%%%%%%%%%%%%

% Determine if you want extra plot for travelling wave of each species
PlotTravelWave = 'y'; % Plot travelling waves
createClusterPlots = 'n';
PlotCommunity = 'n';
CombinedPlot = 'y';   % Plot community

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIM PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create Spacial variables
dx = p.l*4/p.L;                 % define dx based on L
x = linspace(-2*p.l,2*p.l-dx,p.L);    % create continous space
PAD = (abs(x)<=p.l);        % padding outside []

% Space is constructed as
% [2/8 PAD] [1/8 ExtensionBenign] [2/8 Experiment Model] ...
% ... [1/8 ExtensionHostile] [2/8 PAD]

% Indices of when new zone starts
p.PAD1Start = 1;             % Start of left PAD
p.ExtBenStart = (2/8*p.L)+1; % Start Benign Extension
p.ExpModStart = (3/8*p.L)+1; % Start Experimental model
p.ExtHosStart = (5/8*p.L)+1; % Start Hostile Extension
p.PAD2Start =   (6/8*p.L)+1; % Start right PAD

%%%%%%%%%%%%%%%%
%% Simulation %%
%%%%%%%%%%%%%%%%

% Add to Parameter struct
p.PlotTravelWave = PlotTravelWave;
p.CombinedPlot = CombinedPlot;
p.xProj = ((x+p.l/2).*(1/p.l));

%%%%%%%%%%%%%%%%%%%%%%%%
%% Species Generation %%
%%%%%%%%%%%%%%%%%%%%%%%%
p.z = repmat(-p.l/2,p.Species,1); % Range of optimal positions

LastWaves = 5; % How many waves are shown in the plot
p.PlotHandle = zeros(p.clusters,LastWaves,p.L); % Predefine array that plots last (LastWaves) Travelling waves
p.clusterDensity = zeros(p.clusters,p.t_iter);
p.clusterWave = zeros(p.clusters,p.L,p.t_iter);

if length(p.sigma) > 1 % If species specific sigma is used
    p.sig = p.sigma'; % Set each sigma_i as value
else % Otherwise -->  one global sigma
    p.sig = repmat(p.sigma(1),1,p.Species)';
end

p.R = p.R';
p.T = p.T';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get immigration rates %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.I = p.sig;
N = zeros(p.Species, p.L); % Set densities initially to 0

if length(p.alphaIntra) <= p.clusters % Modified competition matrix
    p.comp = repmat(p.alphaInter, p.Species); % Initially all to interspecific
    for i = 1:p.Species
        p.comp(i,i) = p.alphaIntra; % change diagonal values to intraspecific
    end
elseif length(p.alphaIntra) > p.clusters % If custom comp matrix is true
    p.comp = reshape(p.alphaIntra, p.clusters, p.clusters); % Extract into competition matrix
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate growth and dispersal %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.r = 1./(sqrt(2*pi.*p.T.^2)).*exp((-abs(x-p.z).^2)./(2.*p.T)).*(p.T*sqrt(2*pi)).*p.R;
% Extend growth rate beyond experimental scope on benign end
for i = 1:p.clusters
    p.r(i, x < p.z(i)) = p.R(i);
end

%%% Calculate dispersal kernel %%%
k = 1./sqrt(2.*p.sig.^2).* exp(-sqrt(2./(p.sig.^2))*abs(x)); % Laplace Kernel (Lutscher 2019, p.18, Eq. 2.27)
fk = zeros(size(N,1),p.L); % Predefine fft(k)
for s = 1:size(N,1)
    fk(s,:) = fft(k(s,:)); % calculate fft(k)
end

%%%%%%%%%%%%%%%%%%
%% SIMULATE RUN %%
%%%%%%%%%%%%%%%%%%

SpeciesClean = 'n'; % Select if dead species should be deleted

% Calculate original density
clusterNames = ["Spa", "midM", "ArtFes", "Ely"];
if p.clusters == 3
    clusterNames = ["Spa", "midM", "Ely"];
end


%% Time-iterations
for t = 1:p.t_iter
    p.t = t;
    
    p.OriginalDensity = zeros(length(clusterNames), size(p.control.x,2));
    
    cidx = 0;
    for clusterName = clusterNames
        cidx = cidx +1;
        if length(p.betaDens) > 1 % if beta is species specific
            p.OriginalDensity(cidx,:) = -(log(1-p.control.(clusterName)(t,:))/p.betaDens(cidx));
        else % Otherwise, global beta
            p.OriginalDensity(cidx,:) = -(log(1-p.control.(clusterName)(t,:))/p.betaDens);
        end
    end
    
    % Extension with same density as last density in model space
    for i = 1:p.clusters
        p.OriginalDensity(i,1:50) = repmat(p.OriginalDensity(i,51),1,50);
        p.OriginalDensity(i,151:end) = repmat(p.OriginalDensity(i,150),1,50);
    end
    
    % Extension with both ends declining to 0
    % for i = 1:p.clusters
    %     p.OriginalDensity(i,1:50) = linspace(0,p.OriginalDensity(i,51),50);
    %     p.OriginalDensity(i,151:end) = linspace(p.OriginalDensity(i,150),0,50);
    % end
    
    
    %% Interpolate Experimental data to match model resolution
    % Original x-values (uniform spacing)
    x_old = linspace(1, 100, 200); % Original indices (from 1 to 100)
    x_new = linspace(1, 100, p.L*(4/8)); % New x-values (4/8 of model are nonPad space)
    
    % Preallocate the new matrix for interpolated data
    InterpolatedDensity = zeros(size(p.OriginalDensity,1), p.L*(4/8)); % 4/8 of model are nonPad space
    
    % Interpolate each row
    for i = 1:size(p.OriginalDensity,1)
        InterpolatedDensity(i, :) = interp1(x_old, p.OriginalDensity(i, :), x_new, 'linear');
    end
    
    % Add padding to Interpolated Density
    pad1Side = zeros(p.clusters,p.L*2/8); % One side of padding with p.clusters rows
    InterpolatedDensity = [pad1Side, InterpolatedDensity, pad1Side]; % Add padding on both sides
    
    % Update p.OriginalDensity to the interpolated version, flip to match
    % benign left, hostile right side
    p.OriginalDensity = fliplr(InterpolatedDensity);
    
    %%% Calculate Immigration %%%
    % FFT over OriginalDensity --> Only dispersal from outside
    f_im = p.r.*p.OriginalDensity.*exp(-p.comp*p.OriginalDensity);
    N_inv = dx*real( fftshift(ifft(fft(f_im,[],2).*fk,[],2),2));    
    p.gamma = p.gamma';
    N_inv = p.gamma.*p.sig.*N_inv;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% New Density calculation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [N,p] = step(N,x, p, fk, PAD); % Growth+Dispersal
    
    N = N + N_inv;  % Add contribution of invasion
    %%%%%%%%%%%%%%%%%%
    
    for i = 1:p.clusters
        N(i,p.ExtBenStart:p.ExpModStart-1) = repmat(N(i,p.ExpModStart),1,p.L/8);
        N(i,p.ExtHosStart:p.PAD2Start-1) = repmat(N(i,p.ExtHosStart-1),1,p.L/8);
    end
    
    ExtThreshold = 0.01;
    if SpeciesClean == 'y'
        [p,N,Species,fk] = Extinction(p,N,ExtThreshold,p.MethodMeta,fk);
    end
    
    % Plot densities during simulation
    if p.CombinedPlot == 'y'
        
        colorsR = [zeros(1,round(p.Species/2)), linspace(0,1,round(p.Species/2))];
        colorsG = [linspace(0,1,round(p.Species/2)), linspace(1,0,round(p.Species/2))];
        colorsB = [linspace(1,0,round(p.Species/2)), zeros(1,round(p.Species/2))];
        colors = [colorsR; colorsG; colorsB];
        
        if PlotCommunity == 'y'
            figure(1)
            subplot(3,3,1:3)
            
            for i = 1:size(N,1)
                if strcmp(p.MethodMeta, 'selected')
                    plot(x,N(i,:), 'Color', colors(:,i))
                elseif strcmp(p.MethodMeta, 'curve')
                    plot(x,N(i,:), 'Color', colors(:, 100-round( 99 * (p.T(i)-p.Tmin) / (p.Tmax-p.Tmin))))
                end
                hold on
            end
            hold off
        end
        
        % Plot of travelling waves
        if p.PlotTravelWave == 'y'
            
            % Refresh travelling waves
            for i = 1:p.clusters
                p.clusterWave(i,:,t) = N(i,:);
                p.clusterDensity(i,t) = sum(N(i,:));
            end
            
            colors = [linspace(0.8,0.0,size(p.PlotHandle,2)); linspace(0.8,0, size(p.PlotHandle,2)); linspace(0.8,0.0, size(p.PlotHandle,2))];
            colorsRGB = [[linspace(0,0,p.clusters/2),linspace(0,1,p.clusters/2)];...
                [linspace(0,1,p.clusters/2),linspace(0.8,0,p.clusters/2)];...
                [linspace(1,0,p.clusters/2),linspace(0,0,p.clusters/2)]];
            colorsRGB = colorsRGB';
            % Plot travelling wave of all Pioneer SM Species
            subplotIdx = 0;
            if createClusterPlots == 'y'
                for i = 1:p.clusters
                    hold on
                    subplotIdx = subplotIdx +1;
                    subplot(3,p.clusters,subplotIdx)
                    for wave = 1:LastWaves
                        if sum(p.PlotHandle(i,wave,:)) == 0
                            break
                        end
                        plot(p.xProj,squeeze(flipud(p.PlotHandle(i,wave,:))), 'Color', colors(:,end+1-wave));
                    end
                    plot(p.xProj,squeeze(flipud(p.PlotHandle(i,1,:))), 'Color', colors(:,end));
                    title(sprintf('Cluster %i', i))
                    xlabel('x')
                    ylabel('\Sigma N_{i}')
                    xlim([0, 1])
                    ylim([0,2])
                end
                hold off
                
                % Cluster density over time
                for i = 1:p.clusters
                    hold on
                    subplot(3,p.clusters,p.clusters+1:p.clusters*2)
                    plot(1:p.t, p.clusterDensity(i,1:p.t),'Color',colorsRGB(i,:))
                    xlabel('t')
                    ylabel('\Sigma N_{i}')
                end
                hold off
                
                % Plot whole community
                subplot(3,p.clusters,(p.clusters*2)+1:p.clusters*3)
                cla
                hold on
                for i = 1:size(N,1)
                    plot(p.xProj,N(i,:))
                    xlabel('x')
                    ylabel('N_{i}')
                    xlim([0, 1])
                end
                hold off
                
            end
        end
        %drawnow
    end
    
end

% Save cluster density over time + Outside of buffer
p.CDensityTime = sum(p.clusterWave(:,:,1:p.t_iter),2);

end

%%%%%%%%%%%%%%%
%% FUNCTIONS %%
%%%%%%%%%%%%%%%

%%% Main function %%%

function [N,p] = step(N, x, p, fk, PAD)
dx = x(2)-x(1);                 % calculate dx

% Calculate next density

% Growth + Dispersal
f = p.r.*N.*exp(-p.comp*N); % Ricker growth
N_next = dx*real( fftshift(ifft(fft(f,[],2).*fk,[],2),2)); % Dispersal

N_next = N_next.*PAD;% Delete densities in PAD area to avoid periodicity
N = N_next;

end