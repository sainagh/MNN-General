function [coorddeformed,Kdof,Fdof,udof]=FEM_nonlinear(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x)

% author: Pietro Sainaghi
% original program written by Erwin Mueller and Ryan Lee

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
DOFFinal = FEMStruct.DOFFinal; % [] indices of DOF that are not grounded
F = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior
k_base = FEMStruct.k_base; % [6,6] stiffness matrix of individual element
R6 = FEMStruct.R6;
RT6 = FEMStruct.RT6;
DOF = FEMStruct.DOF; % [1] number of degrees of freedom in system
Degrees_per_element = FEMStruct.Degrees_per_element;
Final = FEMStruct.Final; % [NDOFnodes] % indices of nodes that are free to move

% design features structure - LinkPropertiesStruct
nonLinearityType = LinkPropertiesStruct.nonLinearityType;  % [1] link length at rest

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
    % assign imaginary value to nonlinear elements of stiffness matrix
    k(:,:,e)=k_base(:,:,e)*diag([1i*x(e) 1 1 1i*x(e) 1 1]); % USE mtimesx?
    %Combine stiffness matrix for solving
    K(Degrees_per_element(e,:),Degrees_per_element(e,:))=K(Degrees_per_element(e,:),Degrees_per_element(e,:))+RT6(:,:,e)*k(:,:,e)*R6(:,:,e);
end
% take out inversion to use sparse instead
% Kinv=K(Final,Final)^-1; %invert the stiffness Matrix

% initialize displacement vectors
U=(zeros(DOFFinal,Ncases));
U_linear=(zeros(DOFFinal,Ncases));
U_nonlinear=(zeros(DOFFinal,Ncases));
U2=(zeros(Ncoord,DOI,Ncases)); %Translation only interesting
U2_linear=(zeros(Ncoord,DOI,Ncases));
U2_nonlinear=(zeros(Ncoord,DOI,Ncases));

% complex stiffness matrices to be used in computation
Kcomplex = K(Final,Final); % remove grounded elements
Kinvcomplex = inv(Kcomplex); % invert matrix

% handle nonlinear matrix components
K_nonlinear = imag(Kcomplex);
Kinv_nonlinear = imag(Kinvcomplex);
% switch between various nonlinear profiles
switch nonLinearityType

    % quadratic profile for beam stiffness
    case 'Quadaratic'
        for j=1:Ncases
            % inverse quadratic on stiffness matrix 
            Kinv_quadratic = sign(Kinv_nonlinear).*sqrt(abs(Kinv_nonlinear));
            % inverse quadratic on force vector
            Fvec = F(Final,j);
            F_quadratic = sign(Fvec).*sqrt(abs(Fvec));
            % compute displacement
            U_nonlinear(:,j)=Kinv_quadratic*F_quadratic;
            i2=1;
            for i=DOFnodes %loop through the nodes that are not fixed
                U2_nonlinear(i,:,j)=transpose(U_nonlinear(DOI*(i2-1)+1:DOI*(i2-1)+DOI,j));
                i2=i2+1;
            end
        end

end

% handle linear matrix components
K_linear = real(Kcomplex); % take real part of complex stiffness matrix
Kinv_linear = real(Kinvcomplex);
for j=1:Ncases
    U_linear(:,j)=Kinv_linear*F(Final,j);
    i2=1;
    for i=DOFnodes %loop through the nodes that are not fixed
        U2_linear(i,:,j)=transpose(U_linear(DOI*(i2-1)+1:DOI*(i2-1)+DOI,j));
        i2=i2+1;
    end
end

% compute displacement from both linear and nonlinear components
for j=1:Ncases
    coorddeformed(:,:,j)=coord_initial+U2_linear(:,:,j)+U2_nonlinear(:,:,j);
end
U = U_linear+U_nonlinear;
U2 = U2_linear+U2_nonlinear;

% extract stiffness matrix
Kdof = Kcomplex;
% extract input forces in dof format
Fdof = F(Final,:);
% extract final displacements in dof format
udof = U;

end
