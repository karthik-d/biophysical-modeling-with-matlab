function [V, C, vAll, cAll, xbox, ybox] = getVoronoiDiagram(x,y,plotIt)
%% FUNCTION to compute and draw Voronoi diagram given input points
% vInBox: vertex positions of all cells (for cells in box interior)
% cInBox: vertex indices per cell (for cells in box interior)

% get number of points
N = length(x);

% get voronoi diagram information
dt              = delaunayTriangulation(x(:),y(:));
[vAll,cAll]     = voronoiDiagram(dt);

% get total number of vertices
NV = size(vAll,1);

%% Remove voronoi cells that connect to the box boundaries

% initialize indices to include
cBoxInd = true(N,1);
vBoxInd = true(NV,1);

% loop over points from 1 to N
for ii = 1:N
    
    % get indices of vertices connected to cell ii
    vert_ind = cAll{ii};
    
    % check whether to include cells
    if ismember(1,vert_ind)
        cBoxInd(ii,1) = 0;
        vBoxInd(vert_ind,1) = false(length(vert_ind),1);
    else
        vx = vAll(vert_ind,1);
        vy = vAll(vert_ind,2);

        if sum(vx<0) >= 1 ||sum(vx>1) >= 1 || sum(vy<0) >= 1 || sum(vy>1) >= 1
            cBoxInd(ii,1) = 0;
            vBoxInd(vert_ind,1) = false(length(vert_ind),1);
        end
    end
end

% output vertex information for cells in box
C = cAll(cBoxInd);
V = vAll;
xbox = x(cBoxInd);
ybox = y(cBoxInd);

% get number of cells in box
Nbox = sum(cBoxInd);

% if printing, plot points and edges of tessellation
if plotIt
    % create figure object
    figure(1); hold on; box on;
    
    % plot cell centers
    plot(xbox,ybox,'ko','markersize',4,'MarkerFaceColor','k');
    
    % loop over cells
    for jj = 1:Nbox
        % get vertices to plot
        to_plot = C{jj};     
        
        % loop over vertices
        for ii = 1:length(to_plot)
            % wrap vertex checking
            ip1 = ii + 1;
            if ii == length(to_plot)
                ip1 = 1;
            end
            
            % plot vertices
            plot([vAll(to_plot(ii),1),vAll(to_plot(ip1),1)],[vAll(to_plot(ii),2),vAll(to_plot(ip1),2)],'r-');
        end
    end
end


end