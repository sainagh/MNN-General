function [x] = analyticalGradientMethod(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, xInit)

% Algorithm by Jiaji Chen to work with mass-spring model, rewritten for finite truss model by Pietro Sainaghi

% this function uses an analytical gradient descent method to optimize a MNN

%% Explanation of variables

% L [1] loss function, MSE multiplied by a scaling term (works for both unitless MSE and unit MSE) 

% dLdu [Nnodes * 3 * Ncases] derivative of loss function with respect to node displacements (all of them), derivatives of non-output nodes are zero
% dLdu_out [Nnodes * 3 * Ncases] derivative of loss function with respect to node displacements (all of them), derivatives of non-output nodes are zero

% u [Nnodes * 3 * Ncases] displacements of each node (x, y, theta of each node, in order, case by case) [beh1x1 beh1y1 beh1theta1 ... beh1xn beh1yn beh1thetan ... behMx1 behMy1 behMtheta1 ... behMxn behMyn behMthetan] for n nodes and M behaviors
% u_out [NinputANDoutput * 2 * Ncases] displacement of output nodes

% K [Nnodes * 3 * Ncases, Nnodes * 3 * Ncases] stiffness matrix for all behaviors at once

%% input unpacking

% optimizer structure - OptimizerDataStruct
% iterMAX = OptimizerDataStruct.AGDmaxIterations;
iterMAX = 1;
learningRate = 0.01; % Learning Rate
Decay_interval = 50; % Learning rate decay interval
Decay_factor = 5; % Learning rate decay factor
% Batch_num = 1; % Minibatch size UNUSED but shows up in original algorithm
MSEscale = 20; % scaling term for MSE used in algorthm

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases; % [1] number of force behaviors to optimize towards
Target = BehaviorStruct.Target; % [NinputANDoutput, 2, Ncases] desired x and y coordinates for each output node
Elongation_max = BehaviorStruct.Elongation_max; % [1] maximum node displacement in output nodes

% lattice structure - LatticeGeometryStruct
connectivity = LatticeGeometryStruct.connectivity; % [Nbeams,2] connectivity of nodes for each link
Outputline = LatticeGeometryStruct.outputNodes; % [NinputANDoutput] indices of output nodes
NinputANDoutput = LatticeGeometryStruct.NinputANDoutput; % [1] number of input and output nodes in a given lattice
coord_initial = LatticeGeometryStruct.coord_initial; % [Nnodes,3] x,y,z coordinates of each node 
Nnodes = size(coord_initial,1); % [1] number of nodes in the lattice
isOutputVec = zeros(Nnodes,1);
isOutputVec(Outputline) = 1;

% finite truss structure - FEMStruct
forceInputs = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior

%% generate K to x inversion function

[K_contribution] = invertStiffnessMatrix2AxialStiffness(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct);


%% initialize time history storage 

% list of loss function value at each iteration
Loss_hist = []; 

% history of stiffness combinations
Stiffness_hist = [];


%% Optmizer Loop

% assign initial condition on stiffness
x = xInit; % column vector

for iter = 1:iterMAX

    % update learning rate based on decay factor
    if mod(iter, Decay_interval == 0)
        learningRate = learningRate/Decay_factor;
    end

    % store stiffness combination for each iteration
    Stiffness_hist = [Stiffness_hist x];

    % compute loss function (scaled MSE)
    % works with both unitless MSE and unit MSE, but may perform differently based on MSEscale
    MSE = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    L = MSEscale^2/2 * MSE;

    % store loss in vector to plot
    Loss_hist = [Loss_hist L];

    % extract FEM properties used in algorithm
    [coorddeformed,~,stiffnessMatrix]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);

    % loop through each target behavior
    for behIDX = 1:Ncases
        % compute displacements
        displacements(:,:,behIDX) = coorddeformed(:,:,behIDX) - coord_initial; % [Nnodes,3,Ncases]

        % compute difference between target and current deformed coordinates
        u_diff(:,:,behIDX) = coorddeformed(Outputline,1:2,behIDX) - Target(:,:,behIDX); % [NinputANDoutput,2,Ncases]

    end

    % initialize loss function derivative vectors
    dLdu = zeros(Nnodes * 3 * Ncases,1);
    dLdu_out = zeros(NinputANDoutput*2*Ncases,1);
    % initialize displacement vector
    u = zeros(Nnodes * 3 * Ncases,1);
    % initialize force vector
    F = zeros(Nnodes * 3 * Ncases,1);
    % initialize stiffness matrix
    K = zeros(Nnodes * 3 * Ncases,Nnodes * 3 * Ncases);

    % package relevant variables
    for behIDX = 1:Ncases % loop through behaviors

        for outIDX = 1:NinputANDoutput % loop through output nodes
            % assign values to loss function derivative
            dLdu_out( (2*NinputANDoutput*(behIDX-1) + 2*(outIDX-1) + 1) : (2*NinputANDoutput*(behIDX-1) + 2*(outIDX-1) + 2) ) = MSEscale^2 .* [u_diff(outIDX,1,behIDX); u_diff(outIDX,2,behIDX)];
            dLdu( 3*Nnodes*(behIDX-1) +  Outputline(outIDX)*3 - 2 ) = MSEscale^2 * u_diff(outIDX,1,behIDX);
            dLdu( 3*Nnodes*(behIDX-1) +  Outputline(outIDX)*3 - 1 ) = MSEscale^2 * u_diff(outIDX,2,behIDX);
        end

        for nodeIDX = 1:Nnodes % loop through all nodes
            % assign values to displacement vector for all nodes
            u( (3*Nnodes*(behIDX-1) + 3*(nodeIDX-1) + 1) : (3*Nnodes*(behIDX-1) + 3*(nodeIDX-1) + 3) ) = [displacements(nodeIDX,1,behIDX); displacements(nodeIDX,2,behIDX); displacements(nodeIDX,3,behIDX)];
        end

        % assign values to force vector for all nodes
        F( (3*Nnodes*(behIDX-1) + 1) : (3*Nnodes*(behIDX-1) + 3*Nnodes) ) = forceInputs(:,behIDX);

        % assign values to stiffness matrix for all vectors
        K( (3*Nnodes*(behIDX-1) + 1) : (3*Nnodes*(behIDX-1) + 3*Nnodes) , (3*Nnodes*(behIDX-1) + 1) : (3*Nnodes*(behIDX-1) + 3*Nnodes) ) = stiffnessMatrix;

    end

    % NOT NECESSARY assuming derivation is correct (see pptx in this folder)

    % % initialize derivative of force wrt displacement (linear)
    % dFdu = zeros(Nnodes * 3 * Ncases, Nnodes * 3 * Ncases); % second order tensor
    % for IDX_j = 1:size(dFdu,1) % index of force
    %     for IDX_l = 1:size(dFdu,2) % index of derivand displacement
    %             % compute derivative of force wrt displacement (linear)
    %             dFdu(IDX_j,IDX_l) = dFdu(IDX_j,IDX_l) + K(IDX_j,IDX_l);
    % 
    %             % TODO add summation for nonlinear terms here when code gets that
    %     end
    % end
    % 
    % 
    % % initialize derivative of force wrt stiffness
    % dFdK = zeros(Nnodes * 3 * Ncases, Nnodes * 3 * Ncases); % second order tensor (actual object is a diagonal third order tensor)
    % for IDX_j = 1:size(dFdK,1) % index of force
    %     for IDX_n = 1:size(dFdK,3) % index of second dimension of derivant stiffness
    %         % compute derivative of force wrt stiffness
    %         dFdK(IDX_j,IDX_n) = dFdK(IDX_j,IDX_n) + u(IDX_n);
    %     end
    % end

    % initialize derivative of loss function wrt stiffness matrix
    dLdK = zeros(Nnodes * 3 * Ncases, Nnodes * 3 * Ncases); % second order tensor
    for IDX_j = 1:size(dLdK,1)
        for IDX_i = 1:size(dLdK,2) % index of second dimension of derivant stiffness
            for IDX_l = 1:length(dLdu) % index of derivative of loss wrt displacement
                % compute derivative of loss function wrt stiffness
                dLdK(IDX_j,IDX_i) = dLdK(IDX_j,IDX_i) + dLdu(IDX_l) * (-(K(IDX_j,IDX_l))^(-1)) * u(IDX_i);
                A = dLdK(IDX_j,IDX_i)
            end
        end
    end

    % update stiffness combination
    K = K - learningRate * dLdK;

    % obtain value of x from K


    


end











