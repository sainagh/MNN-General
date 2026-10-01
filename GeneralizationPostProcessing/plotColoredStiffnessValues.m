function plotColoredStiffnessValues(LinkPropertiesStruct, LatticeGeometryStruct, x, colormap)

% function by Pietro Sainaghi
% this plots the lattice of a mechanical neural network coloring each link
% based on their stiffness value

%% extract relevant values from structure

% link properties
kLinMin = LinkPropertiesStruct.kLinMin; % minimum controllabel stiffness
kLinMax = LinkPropertiesStruct.kLinMax; % maximum controllable stiffness

% lattice properties
Nbeams = LatticeGeometryStruct.Nbeams; % number of beams
coord_initial = LatticeGeometryStruct.coord_initial; % coordinates of each node [Nnodes DOI]
connectivity = LatticeGeometryStruct.connectivity; % connectivity of each link  [Nbeams 2]

%% load color map

load(colormap);
Ncolors = size(mymap2,1);


%% produce plot

figure
hold on

% plot colored beams
for beamIDX = 1:Nbeams
    plot( [ coord_initial( connectivity(beamIDX,1), 1 ) coord_initial( connectivity(beamIDX,2), 1 ) ],...
        [ coord_initial( connectivity(beamIDX,1), 2 ) coord_initial( connectivity(beamIDX,2), 2 ) ],...
        "Color",mymap2(round( (Ncolors-1) / (kLinMax - kLinMin) * (x(beamIDX) - kLinMin) )+1, :),...
        LineWidth=3)
end

%b plot top and bottom ground
groundT = 0.1;
rectangle('Position',[min(coord_initial(:,1))-groundT min(coord_initial(:,2))-groundT max(coord_initial(:,1))-min(coord_initial(:,1))+groundT*2 groundT],...
    'FaceColor','k','EdgeColor','k')
rectangle('Position',[min(coord_initial(:,1))-groundT max(coord_initial(:,2)) max(coord_initial(:,1))-min(coord_initial(:,1))+groundT*2 groundT],...
    'FaceColor','k','EdgeColor','k')

% plot nodes as dots
for nodeIDX = 1:size(coord_initial,1)
    nodePP = plot( coord_initial(nodeIDX,1), coord_initial(nodeIDX,2),...
        'ok');
    nodePP.MarkerSize=6;
    nodePP.MarkerEdgeColor="black";
    nodePP.MarkerFaceColor="white";
end

axis([min(coord_initial(:,1))-0.1 max(coord_initial(:,1))+0.1 min(coord_initial(:,2))-0.1 max(coord_initial(:,2))+0.1])
axis equal