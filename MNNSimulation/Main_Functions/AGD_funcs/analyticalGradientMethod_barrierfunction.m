function [xOut,terminationMethod] = analyticalGradientMethod_barrierfunction(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, xInit)

% Algorithm by Jiaji Chen to work with mass-spring model, rewritten for finite truss model by Pietro Sainaghi

% this function uses an analytical gradient descent method to optimize a MNN



%% input unpacking

% optimizer structure - OptimizerDataStruct
% iterMAX = OptimizerDataStruct.AGDmaxIterations;

learningRate = 5e-4; % Learning Rate
Decay_interval = 400; % Learning rate decay interval
Decay_factor = 1; % Learning rate decay factor
% Batch_num = 1; % Minibatch size UNUSED but shows up in original algorithm
MSEscale = 1; % scaling term for MSE used in algorthm
muB = 1e-6; % scaling term for barrier function

% termination conditions
iterMAX = 200000; % maximum iterations
MSEChangethreshold = 1e-9; % minimum change in MSE
MSEthreshold = 0.001; % target MSE
averageWindow = 400; % window to look back at plateau intervals
averageWindowFrequency = 5000; % frequency of checks for plateau intervals
averageWindowThreshold = 1e-4; % variation threshold at plateau intervals
chaosThreshold = 10;

% plateau protection
checkPlateauWindow = 200;
plateauVarianceThreshold = 1e-3;
LearningRateMultiplier = 2;


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

% nonlinearities
nonlinearStiffness = LinkPropertiesStruct.nonlinearStiffness;
nonLinearityType = LinkPropertiesStruct.nonLinearityType;


%% Initialize length parameters

% NDOFnodes [] number of nodes that are not grounded
NDOFnodes = length(DOFnodes);
% NDOF [] numebr of degrees of freedom
NDOF = NDOFnodes*DOI;



%% generate K to x inversion function

% initialize stiffness matrix derivative
dKdx_matrix = zeros(NDOF,NDOF,length(xInit));

[stiffnessMatrixDerivative] = computeStiffnessMatrixDerivative(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct);
if size(stiffnessMatrixDerivative) == size(dKdx_matrix)
    dKdx_matrix = stiffnessMatrixDerivative;
else
    error('Issue with computing derivative of stiffness wrt tunable beams')
end

% repackage in matrix form
dKdx = zeros(length(xInit),NDOF^2);
for xIDX = 1:size(dKdx_matrix,3)
    for KIDX = 1:size(dKdx_matrix,1)
        dKdx(xIDX, NDOF*(KIDX-1)+1 : NDOF*(KIDX) ) = dKdx_matrix(KIDX,:,xIDX);
    end
end


%% contruct relationship matrix between DOF and output matrices

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
MSE_hist = [1]; 

% list MSE rate
MSE_rate = [1];

% history of stiffness combinations
Stiffness_hist = [xInit];


%% Optmizer Loop

% assign initial condition on stiffness
x = xInit; % column vector

% iteration counter
iter = 1;
% start change in error
MSEchange = 10;

% change in learning rate flag
MultLearningRateFlag = false;

close all
figure('Name','AGD')

while abs(MSEchange) > MSEChangethreshold

    % update iteration
    iter = iter + 1;

    % store stiffness combination for each iteration
    Stiffness_hist = [Stiffness_hist x];

    % compute barrier function
    B = 0;
    barrO = log( kLinMax - mean([kLinMin,kLinMax]) ) + log( mean([kLinMin,kLinMax]) - kLinMin ); % offset term to always make barrier function positive with minimum at zero
    for xIDX = 1:length(x)
        B = B - muB * ( log( kLinMax - x(xIDX) ) + log( x(xIDX) - kLinMin ) - barrO );
    end

    % compute barrier function derivative
    dBdx = zeros(length(xInit),1);
    for xIDX = 1:length(x)
        dBdx(xIDX) = - ( -1 / (kLinMax - x(xIDX)) + 1 / (x(xIDX) - kLinMin) );
    end

    % compute loss function (scaled MSE)
    % works with both unitless MSE and unit MSE, but may perform differently based on MSEscale
    MSE = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    L = MSEscale^2/2 * MSE + B;
    % disp(['MSE: ', num2str(MSE)])

    % store loss in vector to plot
    MSE_hist = [MSE_hist MSE];

    % compute change in MSE
    % MSEchange = abs(MSE_hist(end-1) - MSE_hist(end));
    % if MSEchange < MSEChangethreshold
    %     disp('Termination: Reached error change threshold')
    % end
    MSEchange = MSE_hist(end-1) - MSE_hist(end);
    if abs(MSEchange) < MSEChangethreshold
        disp('Termination: Reached error change threshold')
        terminationMethod = 'Error Change';
    end
    % store MSE rate in vector to plot
    MSE_rate = [MSE_rate MSEchange];

    % extract FEM properties used in algorithm
    [coorddeformed,Kdof,Fdof,udof]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    F = Fdof;
    u = udof;

    % repackage stiffness matrix into vector
    for KIDX = 1:size(dKdx_matrix,1)
        K(NDOF*(KIDX-1)+1 : NDOF*(KIDX) ,1) = Kdof(KIDX,:);
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

    % compute derivative of cost function with respect to output displacements
    dLdu = zeros(NinputANDoutput*size(Target,2),Ncases);
    for behIDX = 1:Ncases
        for outDOFIDX = 1:(NinputANDoutput*size(Target,2))
            dLdu(outDOFIDX,behIDX) = MSEscale^2 * ( u_out(outDOFIDX,behIDX) - t_out(outDOFIDX,behIDX) );
        end
    end

    % compute derivative of force wrt displacement
    dFdu = zeros(NDOF,NDOF);
    % dFdu = Kdof'; % TEST
    if nonlinearStiffness == false
        dFdu = Kdof;
    else
        switch nonLinearityType
            case 'CB'
                for behIDX = 1:Ncases
                    dFdu_beh(:,:,behIDX) = 2 * Kdof * diag(u(:,behIDX).^2);
                end
                dFdu = sum(dFdu_beh,3);
        end
    end

    % compute derivative of force wrt stiffness
    % dFdK = zeros(NDOF,NDOF^2,Ncases);
    % for behIDX = 1:Ncases
    %     for dofIDX = 1:(NDOF)
    %         for kIDX = 1:(NDOF^2)
    %             % identity of entry is nonzero
    %             rowdivisor = ceil(kIDX/NDOF);
    %             if rowdivisor == dofIDX
    %                 % select specific entry in u
    %                 nonzeroIDX = mod(kIDX,NDOF);
    %                 if nonzeroIDX == 0
    %                     nonzeroIDX = NDOF;
    %                 end
    %                 selectionVec = zeros(1,NDOF);
    %                 selectionVec(nonzeroIDX) = 1;
    %                 % compute nonzero entry
    %                 dFdK(dofIDX,kIDX,behIDX) = u(nonzeroIDX,behIDX);
    %             end
    %         end
    %     end
    % end
    % dFdK_old = dFdK;
    dFdK = zeros(NDOF,NDOF^2,Ncases);
    for behIDX = 1:Ncases
        for dofIDX = 1:(NDOF)
            dFdK(dofIDX, ((dofIDX-1)*NDOF+1) : (dofIDX*NDOF) , behIDX) = u(:,behIDX)';
        end
    end
    % if min(min(min(dFdK_old == dFdK))) == 0
    %     disp('packing error')
    % end


    % compute derivative of displacement with respect to stiffness
    dudK = zeros(NDOF,NDOF^2,Ncases);
    for behIDX = 1:Ncases
        dudK(:,:,behIDX) = dFdu\-dFdK(:,:,behIDX);
    end

    % compute derivative of output displacements with respect to stiffnesssi
    duoutdK = zeros(NinputANDoutput*size(Target,2),NDOF^2,Ncases);
    for outIDX = 1:length(outINDOF_indices)
        duoutdK((outIDX-1)*2+1,:) = dudK(outINDOF_indices(outIDX)*3 - 2,:);
        duoutdK((outIDX-1)*2+2,:) = dudK(outINDOF_indices(outIDX)*3 - 1,:);
    end

    % compute derivative of cost function with respect to stiffness
    dLdK_store = zeros(NinputANDoutput*size(Target,2),NDOF^2,Ncases);
    dLdK = zeros(NDOF^2,1,Ncases);
    for behIDX = 1:Ncases
        dLdK_store(:,:,behIDX) = duoutdK(:,:,behIDX).*dLdu(:,behIDX);
        dLdK(:,:,behIDX) = sum(dLdK_store(:,:,behIDX))';
    end

    % compute derivative of cost function with respect to beam value
    dLdx = zeros(length(xInit),1,Ncases);
    for behIDX = 1:Ncases
        dLdx(:,:,behIDX) = dKdx*dLdK(:,:,behIDX) + muB * dBdx;
    end

    % compute unit scaling term
    largestterm = max(nonzeros(dLdx));
    smallestterm = min(nonzeros(dLdx));
    unitscaler = 1/largestterm;
    

    % compute gradient descent
    for behIDX = 1:Ncases
        % x = x - unitscaler * learningRate * dLdx(:,:,behIDX);
        x = x - unitscaler * learningRate/abs(MSEchange) * dLdx(:,:,behIDX);
    end

    % enforce maximum and minimum stiffness
    % should be redundant but used in debugging
    HitBound(iter+1) = 0;
    for xIDX = 1:length(x)
        if x(xIDX) > kLinMax
            x(xIDX) = kLinMax*0.9;
            disp('Hit an upper bound')
            HitBound(iter+1) = MSE;
        elseif x(xIDX) < kLinMin
            x(xIDX) = kLinMin*0.9;
            disp('Hit a lower bound')
            HitBound(iter+1) = MSE;
        end
    end

    subplot(2,1,1)
    hold off
    plot(MSE_hist,'k-')
    hold on
    plot(HitBound,'r*')
    xlabel('Iteration')
    ylabel('MSE')
    title(['MSE Current: ',num2str(MSE),' Min: ',num2str(min(MSE_hist))])
    set(gca, 'YScale', 'log')
    axis([0 iter MSEthreshold max(MSE_hist)])
    drawnow
    subplot(2,1,2)
    plot(MSE_rate,'k-')
    xlabel('Iteration')
    ylabel('MSE Rate')
    axis([0 iter -chaosThreshold chaosThreshold])
    title(['MSE Rate: ',num2str(MSEchange)])
    % set(gca, 'YScale', 'log')
    drawnow

    if iter > iterMAX
        MSEchange = 0;
        disp('Termination: Maximum Number of Iterations')
        terminationMethod = 'Iteration Count';
    end

    if MSE < MSEthreshold
        MSEchange = 0;
        disp('Termination: Target MSE reached')
        terminationMethod = 'Reached Target MSE';
    end

    % if mod(iter, averageWindowFrequency) == 0
    %     valueWindows = MSE_hist((end-averageWindow):end);
    %     MSEAVwindchange = max([abs(max(valueWindows)-mean(valueWindows)) abs(min(valueWindows)-mean(valueWindows))]);
    %     if MSEAVwindchange < averageWindowThreshold
    %         MSEchange = 0;
    %         disp('Termination: Variation in Iteration Window is less than threshold')
    %     end
    % end

    if (iter > 1000) && (MSE_hist(end) > MSE_hist(end-1)) && (abs(MSE_hist(end) - MSE_hist(end-1)) > chaosThreshold)
        MSEchange = 0;
        disp('Termination: Chaos')
        terminationMethod = 'Chaos';
    end

    if (MSE > 10000)
        MSEchange = 0;
        disp('Termination: Went the wrong way')
        terminationMethod = 'Wrong Way';
    end

    % if (mod(iter,checkPlateauWindow) == 0) && (MultLearningRateFlag == false) && (MSE > 0.05)
    %     if abs(max(MSE_rate( (iter-checkPlateauWindow+1) : checkPlateauWindow))) < plateauVarianceThreshold
    %         learningRate = learningRate * LearningRateMultiplier;
    %         MultLearningRateFlag = true;
    %         disp('Plateau Protection')
    %     end
    % end
    % if (mod(iter,checkPlateauWindow) == (checkPlateauWindow*0.5)) && (MultLearningRateFlag == true)
    %     learningRate = learningRate / LearningRateMultiplier;
    %     MultLearningRateFlag = false;
    % end



end

[minErr minIDX] = min(MSE_hist);
xOut = Stiffness_hist(:,minIDX);

disp('Finished Optimization')









