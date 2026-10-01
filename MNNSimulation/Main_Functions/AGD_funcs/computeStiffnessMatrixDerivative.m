function [dKdx] = computeStiffnessMatrixDerivative(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct)


% author: Pietro Sainaghi
% using architecture written by Erwin Muelder and Ryan Lee

% this function evaluates the matrix terms needed to invert a given resulting MNN stiffness matrix
% to compute axial stiffness values for each beam

%% function inputs

% lattice structure - LatticeGeometryStruct
Nbeams = LatticeGeometryStruct.Nbeams; % [1] number of beams in lattice
coord_initial = LatticeGeometryStruct.coord_initial; % [Ncoord,3] x,y,z coordinates of each node 
DOI = LatticeGeometryStruct.DOI; % [1] dimensions in space = 3
Ncoord = LatticeGeometryStruct.Ncoord;

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases;

% finite truss structure - FEMStruct
DOFnodes = FEMStruct.DOFnodes; % [] indices of nodes that are not grounded
DOFFinal = FEMStruct.DOFFinal;
F = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior
k_base = FEMStruct.k_base; % [6,6] stiffness matrix of individual element
R6 = FEMStruct.R6;
RT6 = FEMStruct.RT6;
DOF = FEMStruct.DOF; % [1] number of degrees of freedom in system
Degrees_per_element = FEMStruct.Degrees_per_element;
Final = FEMStruct.Final;

% optimization problem structure - OptimizerDataStruct
% ForceScaling = OptimizerDataStruct.ForceScaling; % [boolean] yes or no to force scaling

%% variable explanation

% K_contribution [DOF,DOF,Nbeams] for each tunable beam, the imaginary terms show what entries are changed by the tuning, used to contruct a derivative matrix
K_contribution = zeros(DOF,DOF,Nbeams);
% KDOF_contribution [NDOFnodes*DOI,NDOFnodes*DOI,Nbeams] for each tunable beam the imaginary terms show what entries are changed by the tuning, used to contruct a derivative matrix
KDOF_contribution = zeros(length(Final),length(Final),Nbeams);
% dKdx [NDOFnodes*DOI,NDOFnodes*DOI,Nbeams] derivative of stiffness matrix wrt tunable beam stiffness values
dKdx = zeros(length(Final),length(Final),Nbeams);

%% create stiffness contribution matrices


% stiffness matrix for individual beam with control stiffness value
k_beam=(zeros(2*DOI,2*DOI,Nbeams));

% stiffness contribution matrix for full lattice
for beamIDX=1:Nbeams
    % set tuned stiffness as imaginary number to isolate it from non-tunable stiffness values
    k_beam(:,:,beamIDX)=k_base(:,:,beamIDX)*diag([sqrt(-1) 1 1 sqrt(-1) 1 1]); 
    % Compute contribution matrix
    K_contribution(Degrees_per_element(beamIDX,:),Degrees_per_element(beamIDX,:),beamIDX)=RT6(:,:,beamIDX)*k_beam(:,:,beamIDX)*R6(:,:,beamIDX);
end

% stiffness contribution matrix for DOF elements of lattice
for beamIDX=1:Nbeams
    KDOF_contribution(:,:,beamIDX) = K_contribution(Final,Final,beamIDX);
end

%% compute derivatives

for beamIDX = 1:Nbeams
    for KxIDX = 1:size(KDOF_contribution,1)
        for KyIDX = 1:size(KDOF_contribution,2)
            dKdx(KxIDX,KyIDX,beamIDX) = imag( KDOF_contribution( KxIDX,KyIDX,beamIDX ) );
        end
    end
end





end


