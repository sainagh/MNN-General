function [x] = analyticalGradientMethod(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, xInit)

% Algorithm by Jiaji Chen to work with mass-spring model, rewritten for finite truss model by Pietro Sainaghi

% this function uses an analytical gradient descent method to optimize a MNN



%% input unpacking

% optimizer structure - OptimizerDataStruct
% iterMAX = OptimizerDataStruct.AGDmaxIterations;
iterMAX = 2000;
learningRate = 0.01; % Learning Rate
Decay_interval = 100; % Learning rate decay interval
Decay_factor = 5; % Learning rate decay factor
% Batch_num = 1; % Minibatch size UNUSED but shows up in original algorithm
MSEscale = 25; % scaling term for MSE used in algorthm
MSEthreshold = 1e-5;
ChainRuleMultiplier = 1;

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
DOI = LatticeGeometryStruct.DOI; % [1] dimensions in plane, corresponding to x, y, thetaz
Nbeams = LatticeGeometryStruct.Nbeams; % [1] number of tunable beams

% finite truss structure - FEMStruct
F_all = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior
DOFnodes = FEMStruct.DOFnodes; % [] indices of nodes that are not grounded

% beam properties structure - LinkPropertiesStruct
kLinMax = LinkPropertiesStruct.kLinMax;
kLinMin = LinkPropertiesStruct.kLinMin;


%% Explanation of variables

% NDOFnodes [] number of nodes that are not grounded
NDOFnodes = length(DOFnodes);
% NDOF [] numebr of degrees of freedom
NDOF = NDOFnodes*DOI;

% dLdu [NinputANDoutput*2,Ncases] stores derivatives of cost wrt output displacements (x and y) for each behavior
dLdu = zeros(NinputANDoutput*size(Target,2),Ncases);
% F [NDOFnodes*DOI,Ncases] stores force and moment inputs at each node, for each behavior
F = zeros(NDOFnodes*DOI,Ncases);
% u [NDOFnodes*DOI,Ncases] stores displacements (x and y) and rotations (in plane) for each node
u = zeros(NDOFnodes*DOI,Ncases);
% K [NDOFnodes*DOI,NDOFnodes*DOI] stiffness matrix
K = zeros(NDOFnodes*DOI,NDOFnodes*DOI);
% dFdu [NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases] derivative of force wrt displacements, defined as F and u above
dFdu = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases);
% dudF [NDOFnodes*DOI,NDOFnodes*DOI] derivative of displacement wrt force, dudF = -(dFdu)^-1, see derivation
dudF = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases);
% dFdK [NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI] derivative of force wrt stiffness, defined as F and K above; data structure is redundant, [NDOFnodes*DOI,NDOFnodes*DOI] is valid as there's zero terms for all times first and second indices are different, see derivation slides
dFdK = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
% dudK [NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI] derivative of displacement wrt stiffness, defined as u and K above
dudK = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
% duoutdK [NinputANDoutput*2,Ncases,NDOFnodes*DOI,NDOFnodes*DOI] derivative of output displacements wrt stiffness, subset of dudK
duoutdK = zeros(NinputANDoutput*size(Target,2),Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
% dLdK [NDOFnodes*DOI,NDOFnodes*DOI] derivative of cost function wrt stiffness, defined as L and K above
dLdK = zeros(NDOFnodes*DOI,NDOFnodes*DOI);
% dKdx [NDOFnodes*DOI,NDOFnodes*DOI,Nbeams] derivative of stiffness wrt control stiffness values
dKdx = zeros(NDOFnodes*DOI,NDOFnodes*DOI,Nbeams);
% dLdx [Nbeams] derivative of cost function wrt control stiffness values
dLdx = zeros(Nbeams,1);

%% generate K to x inversion function

[stiffnessMatrixDerivative] = computeStiffnessMatrixDerivative(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct);
if size(stiffnessMatrixDerivative) == size(dKdx)
    dKdx = stiffnessMatrixDerivative;
else
    error('Issue with computing derivative of stiffness wrt tunable beams')
end

%% contruct relationship matrix between DOF and output mateices

% identify indices in NDOFnodes format corresponding to output nodes
outINDOF_indices = [];
for outNodeIDX = 1:length(Outputline)
    for dofNodeIDX = 1:length(DOFnodes)
        if Outputline(outNodeIDX) == DOFnodes(dofNodeIDX)
            outINDOF_indices = [outINDOF_indices dofNodeIDX];
        end
    end
end


%% initialize time history storage 

% list of loss function value at each iteration
MSE_hist = []; 

% history of stiffness combinations
Stiffness_hist = [];


%% Optmizer Loop

% assign initial condition on stiffness
x = xInit; % column vector

% iteration counter
iter = 1;
% start change in error
MSEchange = 10;

figure('Name','AGD')

while MSEchange > MSEthreshold

    % update iteration
    iter = iter + 1;

    % update learning rate based on decay factor
    % if mod(iter, Decay_interval == 0)
    %     learningRate = learningRate/Decay_factor;
    % end

    % store stiffness combination for each iteration
    Stiffness_hist = [Stiffness_hist x];

    % compute loss function (scaled MSE)
    % works with both unitless MSE and unit MSE, but may perform differently based on MSEscale
    MSE = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    L = MSEscale^2/2 * MSE;
    % disp(['MSE: ', num2str(MSE)])

    % store loss in vector to plot
    MSE_hist = [MSE_hist MSE];

    % extract FEM properties used in algorithm
    [coorddeformed,Kdof,Fdof,udof]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    K = Kdof;
    F = Fdof;
    u = udof;

    % check if an impossible stiffness matrix is generated
    % good for debugging
    if norm(Fdof-Kdof*udof) > 1e-10
        error('Stiffness matrix is wrong');
    end

    % store all displacements, including grounds
    % there are easier faster ways of doing this, but this is good for debugging as it shows all nodes
    u_all = zeros(Nnodes,DOI,Ncases);
    for behIDX=1:Ncases
        nodeIDX = 1;
        for freeIDX=DOFnodes %loop through the nodes that are not fixed
            u_all(freeIDX,:,behIDX)=transpose(udof(DOI*(nodeIDX-1)+1:DOI*(nodeIDX-1)+DOI,behIDX));
            nodeIDX = nodeIDX + 1;
        end
    end

    % restructure output node displacements and target displacements
    u_out = zeros(NinputANDoutput*size(Target,2),Ncases);
    t_out = zeros(NinputANDoutput*size(Target,2),Ncases);
    for behIDX = 1:Ncases
        for UxIDX = 1:length(Outputline)

            % actual displacements
            u_out((UxIDX-1)*2+1,behIDX) = udof( outINDOF_indices(UxIDX)*3 - 2, behIDX );
            u_out((UxIDX-1)*2+2,behIDX) = udof( outINDOF_indices(UxIDX)*3 - 1, behIDX );

            % target displacements
            t_out((UxIDX-1)*2+1,behIDX) = Target(UxIDX,1,behIDX) - coord_initial(Outputline(UxIDX),1);
            t_out((UxIDX-1)*2+2,behIDX) = Target(UxIDX,2,behIDX) - coord_initial(Outputline(UxIDX),2);

            % debugging
            if u_out((UxIDX-1)*2+1,behIDX) ~= u_all(Outputline(UxIDX),1,behIDX)
                error('Mismatch in displacement vector x')
            end
            if u_out((UxIDX-1)*2+2,behIDX) ~= u_all(Outputline(UxIDX),2,behIDX)
                error('Mismatch in displacement vector x')
            end
        end
    end
   
    % compute dLdu
    dLdu = zeros(NinputANDoutput*size(Target,2),Ncases);
    for behIDX = 1:Ncases
        for outDOFIDX = 1:(NinputANDoutput*size(Target,2))
            dLdu(outDOFIDX,behIDX) = MSEscale^2 * ( u_out(outDOFIDX,behIDX) - t_out(outDOFIDX,behIDX) );
        end
    end

    % % compute dFdu (linear)
    % dFdu = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases);
    % for FdofIDX = 1:NDOF
    %     for FbehIDX = 1:Ncases
    %         for udofIDX = 1:NDOF
    %             for ubehIDX = 1:Ncases
    %                 for dummyIDX = 1:NDOF
    %                     dFdu(FdofIDX,FbehIDX,udofIDX,ubehIDX) = dFdu(FdofIDX,FbehIDX,udofIDX,ubehIDX) + K(FdofIDX,dummyIDX) * (dummyIDX == udofIDX) * (FbehIDX == ubehIDX);
    %                 end
    %             end
    %         end
    %     end
    % end
    % 
    % % compute dFdK
    % dFdK = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
    % for FdofIDX = 1:NDOF
    %     for FbehIDX = 1:Ncases
    %         for KdofIDX1 = 1:NDOF
    %             for KdofIDX2 = 1:NDOF
    %                 for dummyIDX = 1:NDOF
    %                     dFdK(FdofIDX,FbehIDX,KdofIDX1,KdofIDX2) = dFdK(FdofIDX,FbehIDX,KdofIDX1,KdofIDX2) + (FdofIDX == KdofIDX1) * (dummyIDX == KdofIDX2) * u(dummyIDX,FbehIDX);
    %                 end
    %             end
    %         end
    %     end
    % end

    % compute dudK
    % dudK = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
    % for ubehIDX = 1:Ncases
    %     for KdofIDX1 = 1:NDOF
    %         for KdofIDX2 = 1:NDOF
    %             for FdofIDX = 1:NDOF
    %                 for FbehIDX = 1:Ncases
    %                     if (ubehIDX == FbehIDX) && (FdofIDX == KdofIDX1)
    %                         D = -K\u(:,ubehIDX);
    %                     else
    %                         D = zeros(NDOF,1);
    %                     end
    %                     dudK(:,ubehIDX,KdofIDX1,KdofIDX2) = dudK(:,ubehIDX,KdofIDX1,KdofIDX2) + D;
    %                 end
    %             end
    %         end
    %     end
    % end

    % compute dudK
    dudK = zeros(NDOFnodes*DOI,Ncases);
    for ubehIDX = 1:Ncases
        dudK(:,ubehIDX) = -K\u(:,ubehIDX);
    end

    % for testI = 1:size(dudK,3)
    %     for testJ = 1:size(dudK,4)
    % 
    %         if norm(abs(dudK(:,:,testI,testJ) - dudK_reduced(:,:))) > 1e-5
    %             disp('not equivalent')
    %         end
    % 
    %     end
    % end


    % compute duoutdK
    % duoutdK = zeros(NinputANDoutput*size(Target,2),Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
    % for outIDX = 1:length(outINDOF_indices)
    %     duoutdK((outIDX-1)*2+1,:,:,:) = dudK(outINDOF_indices(outIDX)*3 - 2,:,:,:);
    %     duoutdK((outIDX-1)*2+2,:,:,:) = dudK(outINDOF_indices(outIDX)*3 - 1,:,:,:);
    % end

    % compute duoutdK
    duoutdK = zeros(NinputANDoutput*size(Target,2),Ncases);
    for outIDX = 1:length(outINDOF_indices)
        duoutdK((outIDX-1)*2+1,:) = dudK(outINDOF_indices(outIDX)*3 - 2,:);
        duoutdK((outIDX-1)*2+2,:) = dudK(outINDOF_indices(outIDX)*3 - 1,:);
    end
    
    % compute dLdK
    % dLdK = zeros(NDOFnodes*DOI,NDOFnodes*DOI);
    % for KdofIDX1 = 1:NDOF
    %     for KdofIDX2 = 1:NDOF
    %         for udofIDX = 1:(NinputANDoutput*size(Target,2))
    %             for ubehIDX = 1:Ncases
    %                 dLdK(KdofIDX1,KdofIDX2) = dLdK(KdofIDX1,KdofIDX2) + dLdu(udofIDX,ubehIDX) * duoutdK(udofIDX,ubehIDX,KdofIDX1,KdofIDX2);
    %             end
    %         end
    %     end
    % end

    % compute dLdK
    dLdK = zeros(NDOFnodes*DOI,NDOFnodes*DOI);
    for ubehIDX = 1:Ncases
        dLdK = dLdK;
    end

    dLdK = dLdK .* ChainRuleMultiplier;

    % compute dLdx
    dLdx = zeros(Nbeams,1);
    for xIDX = 1:Nbeams
        for KdofIDX1 = 1:NDOF
            for KdofIDX2 = 1:NDOF
                dLdx(xIDX) = dLdx(xIDX) + dLdK(KdofIDX1,KdofIDX2) * dKdx(KdofIDX1,KdofIDX2,xIDX);
            end
        end
    end
    dLdx = dLdx .* ChainRuleMultiplier;

    % compute gradient descent
    x = x - learningRate * dLdx;

    % enforce maximum and minimum stiffness
    for xIDX = 1:length(x)
        if x(xIDX) > kLinMax
            x(xIDX) = kLinMax;
        elseif x(xIDX) < kLinMin
            x(xIDX) = kLinMin;
        end
    end

    
    plot(MSE_hist)
    xlabel('Iteration')
    ylabel('MSE')
    drawnow

    if iter > iterMAX
        MSEchange = 0;
    end

end

disp('Finished Optimization')









