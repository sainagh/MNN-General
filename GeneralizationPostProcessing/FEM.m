function [coorddeformed,Kf2u,Ku2f]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x)

% author: Pietro Sainaghi
% original program written by Reinier Kuppens and Ryan Lee

% this function evaluates the deformed configuration of the lattice for
% each behavior using the finite truss method

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

%% function outputs

coorddeformed=zeros(Ncoord,DOI,Ncases);


%% computation
Ncoord = size(coord_initial,1);
K=(zeros(DOF));
k=(zeros(2*DOI,2*DOI,Nbeams));
Ncases = size(F,2);
for e=1:Nbeams
    %scale the stiffnes accoding to x
    k(:,:,e)=k_base(:,:,e)*diag([x(e) 1 1 x(e) 1 1]); % USE mtimesx?
    %Combine stiffness matrix for solving
    K(Degrees_per_element(e,:),Degrees_per_element(e,:))=K(Degrees_per_element(e,:),Degrees_per_element(e,:))+RT6(:,:,e)*k(:,:,e)*R6(:,:,e);
end
Kinv=K(Final,Final)^-1; %invert the stiffness Matrix
U=(zeros(DOFFinal,Ncases));
U2=(zeros(Ncoord,DOI,Ncases)); %Translation only interesting

for j=1:Ncases
    U(:,j)=Kinv*F(Final,j);
    i2=1;
    for i=DOFnodes %loop through the nodes that are not fixed
        U2(i,:,j)=transpose(U(DOI*(i2-1)+1:DOI*(i2-1)+DOI,j));
        i2=i2+1;
    end
    
    coorddeformed(:,:,j)=coord_initial+U2(:,:,j);

end
Kf2u = Kinv;
Ku2f = K;
end
