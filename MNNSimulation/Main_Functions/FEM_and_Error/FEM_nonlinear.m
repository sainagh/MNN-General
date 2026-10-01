function [coorddeformed,Kdof,Fdof,udof]=FEM_nonlinear(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x)

% author: Pietro Sainaghi
% original program written by Erwin Muller and Ryan Lee

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
MaxForce = BehaviorStruct.MaxForce;


% link properties structure
nonLinearityType = LinkPropertiesStruct.nonLinearityType;
MaxLinkU = max([LinkPropertiesStruct.MaxLinkElongation, LinkPropertiesStruct.MaxLinkCompression]);

% optimization problem structure - OptimizerDataStruct
% ForceScaling = OptimizerDataStruct.ForceScaling; % [boolean] yes or no to force scaling

%% function outputs

coorddeformed=zeros(Ncoord,DOI,Ncases);


%% computate stiffness matrices

Ncoord = size(coord_initial,1);
K_complex=(zeros(DOF));
k=(zeros(2*DOI,2*DOI,Nbeams));
Ncases = size(F,2);
for e=1:Nbeams
    %scale the stiffnes accoding to x
    k(:,:,e)=k_base(:,:,e)*diag([1i*x(e) 1 1 1i*x(e) 1 1]); % USE mtimesx?
    %Combine stiffness matrix for solving
    K_complex(Degrees_per_element(e,:),Degrees_per_element(e,:))=K_complex(Degrees_per_element(e,:),Degrees_per_element(e,:))+RT6(:,:,e)*k(:,:,e)*R6(:,:,e);
end
% take out inversion to use sparse instead
% Kinv=K(Final,Final)^-1; %invert the stiffness Matrix
U_ic=(zeros(DOFFinal,Ncases));
U=(zeros(DOFFinal,Ncases));
U2_ic=(zeros(Ncoord,DOI,Ncases)); %Translation only interesting
U2=(zeros(Ncoord,DOI,Ncases));

% remove ground points
Ksub_complex = K_complex(Final,Final);

% separate into linear and nonlinear
Ksub_linear = real(Ksub_complex);
Ksub_nonlinear = imag(Ksub_complex);

% linearized stiffness matrix for initial conditions
Ksub_IC = Ksub_linear + Ksub_nonlinear;



%% numerical solver to compute nonlinear system

% loop through behaviors
for behNum=1:Ncases

    % get displacemnt initial condition for given behavior
    U_ic_inter = U_ic(:,behNum);

    % get forces stored as single vector
    F_forresidual = F(Final,behNum);

    % switch for nonlinearity types
    switch nonLinearityType
        case 'ATAN' % F = Klin * u + Knonlin * (2*Fmax/pi) * atan(u)
            % define nonlinear function 
            fun = @(U) nonlinearstiffness_atan(F_forresidual, Ksub_linear, Ksub_nonlinear, MaxForce, MaxLinkU, U);
        case 'RAMP' % F = Klin * u + Knonlin * u * (u>=0)
            % define nonlinear function
            fun = @(U) nonlinearstiffness_ramp(F_forresidual, Ksub_linear, Ksub_nonlinear, U);
        case 'EXP' % F = Klin * u + Knonlin * Fmax * (exp(u)-1)
            % define nonlienar function
            fun = @(U) nonlinearstiffness_exp(F_forresidual, Ksub_linear, Ksub_nonlinear, MaxForce, U);
        case 'CUBE' % F = Klin * u + Knonlin * Fmax * u^3
            % define nonlienar function
            fun = @(U) nonlinearstiffness_cube(F_forresidual, Ksub_linear, Ksub_nonlinear, MaxForce, MaxLinkU, U);
    end

    % nonlinear solver options
    options = optimoptions('fsolve',...
        'Algorithm','trust-region-dogleg',...
        'Display','none');

    % solve nonlinear system
    U_computed = fsolve(fun, U_ic_inter,options);

    % assign value to nonlinear displacement for given behavior
    U(:,behNum) = U_computed;
end


%% compute deformed lattice

% loop through behaviors
for behNum=1:Ncases
    DOFnum2=1;
    for DOFnum=DOFnodes %loop through the nodes that are not fixed
        U2(DOFnum,:,behNum)=transpose(U(DOI*(DOFnum2-1)+1:DOI*(DOFnum2-1)+DOI,behNum));
        DOFnum2=DOFnum2+1;
    end
    coorddeformed(:,:,behNum)=coord_initial+U2(:,:,behNum);
end

% extract stiffness matrix
Kdof = Ksub_complex;
% extract input forces in dof format
Fdof = F(Final,:);
% extract final displacements in dof format
udof = U;

end
