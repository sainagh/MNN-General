function LatticeGeometryStruct = packageLattice(LinkPropertiesStruct, LatticeGeometryStruct)

% this function runs the function to create the lattice mesh matrices, and
% combines the outputs in a convenient structure

if LatticeGeometryStruct.latticeType == 1
    
    % generate triangular lattice
    [Ncoord, coord_initial, connectivity, bound,inputNodes,outputNodes, Nbeams] ...
        = triangleLatticeFunction(LinkPropertiesStruct, LatticeGeometryStruct);
    
elseif LatticeGeometryStruct.latticeType == 2
    % generate square lattice
    % not yet implemented
    [Ncoord, coord_initial, connectivity, bound,inputNodes, outputNodes]...
        = squareLatticeFunction(LinkPropertiesStruct, LatticeGeometryStruct);
else
    % error message
    latticeType = input( 'lattice Type not recognized. Please enter a lattice type, 1 (triangular) 2(square)');
end

% store lattice mesh into structure LatticeGeometryStruct
% number of nodes in lattice
LatticeGeometryStruct.Ncoord = Ncoord;
% number of beams in the lattice
LatticeGeometryStruct.Nbeams = Nbeams;
% coordinates of each node when lattice is undisturbed
LatticeGeometryStruct.coord_initial = coord_initial;
% list of nodes each link connects
LatticeGeometryStruct.connectivity = connectivity;
% indices of ground nodes
LatticeGeometryStruct.bound = bound;
% indices of input nodes
LatticeGeometryStruct.inputNodes = inputNodes;
% indices of output nodes
LatticeGeometryStruct.outputNodes = outputNodes;