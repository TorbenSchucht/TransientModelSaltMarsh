% Simulates occupation of grid
% to estimate conversion of density to frequency

clear all
rng(100)

p = struct(); % Parameter struct

% Domain size
p.xSize = 9;
p.ySize = 9;
p.densityScale = 5;
% Grid-Dimensions
p.xGrid = 9;
p.yGrid = 9;
p.iter = 500; % Number of individuals that are generated
p.plotOption = 'y';

% Struct for species Parameters
s = struct(); 

% Predefine result matrices
indCount = zeros(p.xGrid, p.yGrid);
freqCount = zeros(p.xGrid, p.yGrid);
densCount = zeros(p.xGrid, p.yGrid);
gridCellNr = (p.xGrid*p.yGrid);

totalDensity = zeros(1,p.iter);
gridFrequency = zeros(1,p.iter);
xVect = zeros(1,p.iter);
yVect = zeros(1,p.iter);

% Generate grid lines
x_edges = linspace(0, p.xSize, p.xGrid + 1); % Vertical lines
y_edges = linspace(0, p.ySize, p.yGrid + 1); % Horizontal lines


for i = 1:p.iter % For each generated individual
    s.x = rand()*p.xSize;  % assign x,y --> Coordinates
    s.y = rand()*p.ySize;
    s.N = rand()*p.densityScale; % assign N --> Density
    xVect(i) = s.x;
    yVect(i) = s.y;
    [CellX, CellY] = FindCell(s,p); % Find grid cell, that belongs to x,y coordinates
    indCount(CellX, CellY) = indCount(CellX, CellY) +1; % Add to individual count
    densCount(CellX, CellY) = densCount(CellX, CellY) + s.N; % Add to density
    if freqCount(CellX, CellY) == 0
        freqCount(CellX, CellY) = 1; % If gridcell is not yet counted, set it to 1
    end
    
    
    totalDensity(i) = sum(sum(densCount));
    gridFrequency(i) = sum(sum(freqCount));
    
if p.plotOption == 'y'
    % --- Scatter Plot ---
    subplot(1,3,1)
    scatter(xVect, yVect, '.')
    title('Individual Positions')
    
    % Plot vertical grid lines
    for k = 1:length(x_edges)
        line([x_edges(k) x_edges(k)], [0 p.ySize], 'Color', 'k', 'LineStyle', '--');
    end
    
    % Plot horizontal grid lines
    for k = 1:length(y_edges)
        line([0 p.xSize], [y_edges(k) y_edges(k)], 'Color', 'k', 'LineStyle', '--');
    end
    
    xlim([0 p.xSize])
    ylim([0 p.ySize])
    axis square
    xlabel('x', 'FontSize', 12)
    ylabel('y', 'FontSize', 12)
    
    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
    
    % --- Density Plot ---
%     subplot(1,3,2)
%     imagesc(x_centers, y_centers, densCount');    set(gca,'YDir','normal');
%     %colorbar
%     title('Density')
%         % Plot vertical grid lines
%     for k = 1:length(x_edges)
%         line([x_edges(k) x_edges(k)], [0 p.ySize], 'Color', 'k', 'LineStyle', '--');
%     end 
%     % Plot horizontal grid lines
%     for k = 1:length(y_edges)
%         line([0 p.xSize], [y_edges(k) y_edges(k)], 'Color', 'k', 'LineStyle', '--');
%     end
%     xlim([0 p.xSize])
%     ylim([0 p.ySize])
%     axis square

    % --- Frequency Plot ---
    subplot(1,3,2)
    imagesc(x_centers, y_centers, freqCount'); 
    set(gca,'YDir','normal');
    title('Occupancy Matrix')
     % Plot vertical grid lines
    for k = 1:length(x_edges)
        line([x_edges(k) x_edges(k)], [0 p.ySize], 'Color', 'k', 'LineStyle', '--');
    end 
    % Plot horizontal grid lines
    for k = 1:length(y_edges)
        line([0 p.xSize], [y_edges(k) y_edges(k)], 'Color', 'k', 'LineStyle', '--');
    end
    xlim([0 p.xSize])
    ylim([0 p.ySize])
    axis square
    xlabel('x', 'FontSize' ,12)

    drawnow
end

    
end

subplot(1,3,3)
plot(totalDensity, gridFrequency./gridCellNr);
xlabel('Total Density')
ylabel('Frequency (Grid Cells)')

% Find Function fit
x = totalDensity;
y = gridFrequency./gridCellNr;
alpha0 = 0.1;
model = @(alpha, x) 1-exp(-alpha * x);
alpha_fit = lsqcurvefit(model, alpha0, x, y);
x_fit = linspace(min(x), max(x), 100); % x-values for function plot
y_fit = model(alpha_fit, x_fit); % y-values for function plot
figure(1);
plot(x, y, 'r.', 'DisplayName', 'Data'); hold on; % Artificial data plot
plot(x_fit, y_fit, 'b-', 'DisplayName', 'Fitted Curve'); % Function plot
xlabel('n'); ylabel('f');
legend;
axis square
grid on;
AddLetters2Plots(gcf,{'a)','b)','c)'}, 'VShift', 0.22, 'HShift', -0.05', 'FontSize', 14)


function [CellX,CellY] = FindCell(s,p)
deltax = p.xSize/p.xGrid;
deltay = p.ySize/p.yGrid;
for i = 1:p.xGrid
    for j = 1:p.yGrid
        minx = deltax*i-deltax;
        maxx = deltax*i;
        miny = deltay*j-deltay;
        maxy = deltay*j;
        if s.x > minx && s.x < maxx && s.y > miny && s.y < maxy % If x and y are within border
            CellX = i;
            CellY = j;
            break
        end      
    end
end
end















