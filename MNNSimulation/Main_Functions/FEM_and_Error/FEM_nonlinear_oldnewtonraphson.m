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
Fmax = max(max(F));

% link properties structure
nonLinearityType = LinkPropertiesStruct.nonLinearityType;
MaxLinkElongation = LinkPropertiesStruct.MaxLinkElongation;

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
[U,S,V] = svd(Ksub_nonlinear)

% linearized stiffness matrix for initial conditions
Ksub_IC = Ksub_linear + Ksub_nonlinear;

%% compute initial condition for newton raphson

%treat system as fully linear

% loop through behaviors
% for behNum=1:Ncases
%     % sparse matrix inversion
%     U_ic(:,behNum)=Ksub_IC\F(Final,behNum);
%     % DOFnum2=1;
%     % for DOFnum=DOFnodes %loop through the nodes that are not fixed
%     %     U2_ic(DOFnum,:,behNum)=transpose(U_ic(DOI*(DOFnum2-1)+1:DOI*(DOFnum2-1)+DOI,behNum));
%     %     DOFnum2=DOFnum2+1;
%     % end
%     %
%     % coorddeformed(:,:,behNum)=coord_initial+U2_ic(:,:,behNum);
% 
% end

%% newton raphson to compute nonlinear system

% loop through behaviors
for behNum=1:Ncases

    % switch based on nonlinear function
    switch nonLinearityType
        case 'CB'
            % get displacemnt initial condition for given behavior
            Uiter = U_ic(:,behNum);
            % start flag for convergence
            deltaNorm = 1000;
            iterNR = 0;
            while deltaNorm > 1e-10
                iterNR = iterNR+1;
                if iterNR > 100
                    error('Newton-Raphson did not converge')
                end
                % factor for input multiplier
                inputmult = 0.0025*2/pi;
                % compute newton raphson function
                L = -F(Final,behNum) + Ksub_linear * Uiter + inputmult .* Ksub_nonlinear * (Uiter).^3;
                % compute newton raphson jacobian
                dLdU = Ksub_linear + 3 * inputmult .* Ksub_nonlinear * diag( (Uiter).^2 );
                % compute updated iteration
                delta = dLdU\L;
                deltaNorm = norm(delta);
                Uiter = Uiter - delta;
            end
            % assign value to nonlinear displacement for given behavior
            U(:,behNum) = Uiter;

        case 'ATAN'
            % get displacemnt initial condition for given behavior
            Uiter = U_ic(:,behNum);
            % start flag for convergence
            deltaNorm = 1000;
            Lnorm = 1000;
            iterNR = 0;
            while ( deltaNorm > 1e-10 ) || ( Lnorm > 1e-10 )
                iterNR = iterNR+1;
                if iterNR > 100
                    error('Newton-Raphson did not converge')
                end
                % factor for input multiplier
                inputmult =  1*2/pi;
                % compute newton raphson function
                L = -F(Final,behNum) + Ksub_linear * Uiter + inputmult .* Ksub_nonlinear * atan(Uiter);
                % compute newton raphson jacobian
                dLdU = Ksub_linear + 3 * inputmult .* Ksub_nonlinear * diag( 1 / ( 1 + Uiter.^2) );
                % compute updated iteration
                delta = dLdU\L;
                deltaNorm = norm(delta);
                Lnorm = norm(L);
                Uiter = Uiter - delta;
            end
            % assign value to nonlinear displacement for given behavior
            U(:,behNum) = Uiter;

        case 'ERF'
            % get displacemnt initial condition for given behavior
            Uiter = U_ic(:,behNum);
            % start flag for convergence
            deltaNorm = 1000;
            iterNR = 0;
            while deltaNorm > 1e-10
                iterNR = iterNR+1;
                if iterNR > 100
                    error('Newton-Raphson did not converge')
                end
                % compute newton raphson function
                L = -F(Final,behNum) + Ksub_linear * Uiter + Ksub_nonlinear * ( erf( Uiter ) );
                % compute newton raphson jacobian
                dLdU = Ksub_linear + Ksub_nonlinear * diag( (1-1e-10) * exp( -( Uiter ).^2 ) );
                % compute updated iteration
                delta = dLdU\L;
                deltaNorm = norm(L);
                Uiter = Uiter - delta;
            end
            % assign value to nonlinear displacement for given behavior
            U(:,behNum) = Uiter;
    end
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
