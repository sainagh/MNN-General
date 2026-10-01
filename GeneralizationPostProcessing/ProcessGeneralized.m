clear
clc
close all

%% load data

[fileList, dataPath] = uigetfile('*.mat','Select Result Structure','MultiSelect', 'on');
load([dataPath,fileList]);

load('highAmp_densesweep.mat')

%% Plot Target Behaviors

% number of output nodes
Noutput = LatticeGeometryStruct.NinputANDoutput;

% discrete color map
TargetBehColorMap = [linspace(1,0,BehaviorStruct.NcasesArray); linspace(0,1,BehaviorStruct.NcasesArray); zeros(1,BehaviorStruct.NcasesArray)]';

figure
hold on
for outIDX = 1:Noutput
        subplot(Noutput,1, Noutput+1-outIDX )
        hold on
        xTAR(:,1) = BehaviorStruct.Target(outIDX,1,:);
        yTAR(:,1) = BehaviorStruct.Target(outIDX,2,:);
        for behIDX = 1:size(xTAR,1)
            plot(xTAR(behIDX,1),yTAR(behIDX,1),'*','Color',TargetBehColorMap(behIDX,:),'LineWidth',1.5,'MarkerSize',9)
        end
        axis([mean(BehaviorStruct.Target(outIDX,1,:))-0.004 mean(BehaviorStruct.Target(outIDX,1,:))+0.004 BehaviorStruct.Target(outIDX,2,1)-0.0005 BehaviorStruct.Target(outIDX,2,1)+0.0005])
        box on
        set(gca,'lineWidth',2.5)
        set(gca,'YTick',[BehaviorStruct.Target(outIDX,2,1)])
        set(gca,'XTick',[mean(BehaviorStruct.Target(outIDX,1,:))-0.004 mean(BehaviorStruct.Target(outIDX,1,:))-0.002 mean(BehaviorStruct.Target(outIDX,1,:)) mean(BehaviorStruct.Target(outIDX,1,:))+0.002 mean(BehaviorStruct.Target(outIDX,1,:))+0.004])
end

%% Plot Continuous Target Space

% get sinusoid amplitude for general behavior
sinAmp = DenseSweepStruct.dx;

% continuous color map
ContinuousBehColorMap = [linspace(1,0,DenseSweepStruct.Ncases); linspace(0,1,DenseSweepStruct.Ncases); zeros(1,DenseSweepStruct.Ncases)]';

% input force angles
FSinSPace = linspace(0,pi,DenseSweepStruct.Ncases);

% discrete input force angles
FSinSPaceDisc = linspace(0,pi,size(xTAR,1));

figure
hold on
% plot
for FIDX = 1:length(FSinSPace)
    plot(FSinSPace(FIDX),1000*sinAmp*cos(FSinSPace(FIDX)),'.','Color',ContinuousBehColorMap(FIDX,:))
end

for FIDX = 1:length(FSinSPaceDisc)
    plot(FSinSPaceDisc(FIDX),1000*sinAmp*cos(FSinSPaceDisc(FIDX)),'*','Color',TargetBehColorMap(FIDX,:),'MarkerSize',15,'LineWidth',2.5)
end

box on
set(gca,'lineWidth',2.5)
axis([0 pi -1.5*sinAmp*1000 1.5*sinAmp*1000])
set(gca,'FontSize',18)

%% Compute Discrete Points MSE and Displacements

% final stiffness combination from study
xFinal = (SweepOutput{1,1}.x)';

% set number of behaviors to 1 for each via point, as MSE is computed
% individually for each
DiscretePointsBehaviorStruct.Ncases = 1;

% max allowed elongation from link structure
DiscretePointsBehaviorStruct.Elongation_max = LinkPropertiesStruct.MaxLinkElongation;

% number of discrete points
Ndiscrete = (BehaviorStruct.NcasesArray);

% compute MSE for each discrete point
for discreteIDX = 1:Ndiscrete

    % set single target for MSE computation
    DiscretePointsBehaviorStruct.Target = BehaviorStruct.Target(:,:,discreteIDX);

    % set single force for MSE computation
    DiscretePointsBehaviorStruct.forces = BehaviorStruct.forces(:,:,discreteIDX);

    % package FEM
    [DOFFinal,Final, F, DOFnodes,k_base,R6,RT6,DOF,Degrees_per_element] = FEM_setup(LinkPropertiesStruct,LatticeGeometryStruct,DiscretePointsBehaviorStruct);
    FEMStructDisc.DOFFinal=DOFFinal;
    FEMStructDisc.Final=Final;
    FEMStructDisc.F=F;
    FEMStructDisc.DOFnodes=DOFnodes;
    FEMStructDisc.k_base=k_base;
    FEMStructDisc.R6=R6;
    FEMStructDisc.RT6=RT6;
    FEMStructDisc.DOF=DOF;
    FEMStructDisc.Degrees_per_element=Degrees_per_element;

    % compute MSE for discrete point
    errorDiscrete(discreteIDX) = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, DiscretePointsBehaviorStruct,FEMStructDisc,OptimizerDataStruct,xFinal);
    
    % compute node displacements for discrete point
    deformedDiscrete(:,:,:,discreteIDX)=FEM(LinkPropertiesStruct, LatticeGeometryStruct, DiscretePointsBehaviorStruct,FEMStructDisc,OptimizerDataStruct,xFinal);
end

%% Compute 'Continuous' Space MSE and Displacements

% set number of behaviors to 1 for each via point, as MSE is computed
% individually for each
ContinuousBehaviorStruct.Ncases = 1;

% max allowed elongation from link structure
ContinuousBehaviorStruct.Elongation_max = LinkPropertiesStruct.MaxLinkElongation;

% number of 'continuous' points
Ncontinuous = (DenseSweepStruct.Ncases);

% compute MSE for each continuous point
for continuousIDX = 1:Ncontinuous

    disp(['ErrDispCont Iter: ',num2str(continuousIDX)])
    
    % set single target for MSE computation
    ContinuousBehaviorStruct.Target = DenseSweepStruct.Target(:,:,continuousIDX);

    % set single force for MSE computation
    ContinuousBehaviorStruct.forces = DenseSweepStruct.forces(:,:,continuousIDX);

    % package FEM
    [DOFFinal,Final, F, DOFnodes,k_base,R6,RT6,DOF,Degrees_per_element] = FEM_setup(LinkPropertiesStruct,LatticeGeometryStruct,ContinuousBehaviorStruct);
    FEMStructVia.DOFFinal=DOFFinal;
    FEMStructVia.Final=Final;
    FEMStructVia.F=F;
    FEMStructVia.DOFnodes=DOFnodes;
    FEMStructVia.k_base=k_base;
    FEMStructVia.R6=R6;
    FEMStructVia.RT6=RT6;
    FEMStructVia.DOF=DOF;
    FEMStructVia.Degrees_per_element=Degrees_per_element;

    % compute MSE for continuous point
    errorContinuous(continuousIDX) = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, ContinuousBehaviorStruct,FEMStructVia,OptimizerDataStruct,xFinal);
     % compute node displacements for discrete point
    deformedContinuous(:,:,:,continuousIDX)=FEM(LinkPropertiesStruct, LatticeGeometryStruct, ContinuousBehaviorStruct,FEMStructVia,OptimizerDataStruct,xFinal);
end

%% plot homogeneous lattice

plotColoredStiffnessValues(LinkPropertiesStruct, LatticeGeometryStruct, zeros(size(xFinal,1), size(xFinal,2)), 'colormap.mat')
axis off

%% plot color map

plotColoredStiffnessValues(LinkPropertiesStruct, LatticeGeometryStruct, xFinal, 'colormap.mat')
axis off

%% plot error curve

figure
hold on
thetaDisc = linspace(0,pi,Ndiscrete);
for discIDX = 1:length(thetaDisc)
    plot(thetaDisc(discIDX),errorDiscrete(discIDX),'d','MarkerSize',8, 'Color',TargetBehColorMap(discIDX,:),'LineWidth',2.5)
end
% plot(thetaDisc, mean(errorDiscrete)*ones(1,Ndiscrete),'d','MarkerSize',8, 'Color',TargetBehColorMap(round(BehaviorStruct.NcasesArray/2),:),'LineWidth',2.5)
thetaCont = linspace(0,pi,Ncontinuous);
% for contIDX = 1:length(thetaCont)
%     plot(thetaCont(contIDX),errorContinuous(contIDX),'.', 'Color',ContinuousBehColorMap(contIDX,:))
% end
% plot(thetaCont,mean(errorContinuous)*ones(1,Ncontinuous),'.', 'Color',ContinuousBehColorMap(round(DenseSweepStruct.Ncases/2),:))
box on
set(gca,'lineWidth',2.5)
axis([0 pi 0.9*min(errorContinuous) 1.1*max(errorContinuous)])
set(gca,'FontSize',18)

mean(errorDiscrete)
mean(errorContinuous)

%% plot deformed coordinates

% extract coordinates of output nodes
OutDeformedDiscrete = deformedDiscrete(LatticeGeometryStruct.outputNodes,:,1,:);
OutDeformedContinuous = deformedContinuous(LatticeGeometryStruct.outputNodes,:,1,:);

% initial coordinates
OutCoordInitial = LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes,:);

% plot all output node displacement in same figure
figure
hold on
for outIDX = 1:Noutput
    % compute displacements
    for discreteIDX = 1:Ndiscrete
        OutDispDiscrete(outIDX,:,1,discreteIDX) = OutDeformedDiscrete(outIDX,:,1,discreteIDX)-OutCoordInitial(outIDX,:);
    end
    for continuousIDX = 1:Ncontinuous
        OutDispContinuous(outIDX,:,1,continuousIDX) = OutDeformedContinuous(outIDX,:,1,continuousIDX)-OutCoordInitial(outIDX,:);
    end

    % plot displacements
    XdispD(:) = OutDispDiscrete(outIDX,1,1,:);
    YdispD(:) = OutDispDiscrete(outIDX,2,1,:);
    plot(XdispD,YdispD,'*b')
    plot(-XdispD,-YdispD,'*b')
    XdispC(:) = OutDispContinuous(outIDX,1,1,:);
    YdispC(:) = OutDispContinuous(outIDX,2,1,:);
    plot(XdispC,YdispC,'-b')
    plot(-XdispC,-YdispC,'-b')
end

% plot all output node displacement in subplots
figure
hold on
for outIDX = 1:Noutput
        subplot(Noutput,1, Noutput+1-outIDX )
        hold on
        xTAR(:,1) = BehaviorStruct.Target(outIDX,1,:);
        yTAR(:,1) = BehaviorStruct.Target(outIDX,2,:);
        for behIDX = 1:size(xTAR,1)
            plot(xTAR(behIDX,1),yTAR(behIDX,1),'*','Color',TargetBehColorMap(behIDX,:),'LineWidth',1.5,'MarkerSize',9)
        end
        xDIS(:,1) = OutDeformedDiscrete(outIDX,1,1,:);
        yDIS(:,1) = OutDeformedDiscrete(outIDX,2,1,:);
        for behIDX = 1:size(xTAR,1)
            plot(xDIS(behIDX,1),yDIS(behIDX,1),'d','Color',TargetBehColorMap(behIDX,:),'LineWidth',1.5,'MarkerSize',9)
        end
        xCON(:,1) = OutDeformedContinuous(outIDX,1,1,:);
        yCON(:,1) = OutDeformedContinuous(outIDX,2,1,:);
        for behIDX = 1:size(xCON,1)
            plot(xCON(behIDX,1),yCON(behIDX,1),'.','Color',ContinuousBehColorMap(behIDX,:))
        end
        axis([mean(BehaviorStruct.Target(outIDX,1,:))-0.002 mean(BehaviorStruct.Target(outIDX,1,:))+0.002 BehaviorStruct.Target(outIDX,2,1)-0.0005 BehaviorStruct.Target(outIDX,2,1)+0.0005])
        box on
        set(gca,'lineWidth',2.5)
        set(gca,'YTick',[BehaviorStruct.Target(outIDX,2,1)])
        set(gca,'XTick',[mean(BehaviorStruct.Target(outIDX,1,:))-0.004 mean(BehaviorStruct.Target(outIDX,1,:))-0.002 mean(BehaviorStruct.Target(outIDX,1,:)) mean(BehaviorStruct.Target(outIDX,1,:))+0.002 mean(BehaviorStruct.Target(outIDX,1,:))+0.004])

end

%% Compute Reduced Transformation Matrix

% set number of behaviors to 1 for each decoupled point, as MSE is computed
% individually for each
DecoupledBehaviorStruct.Ncases = 1;

% max allowed elongation from link structure
DecoupledBehaviorStruct.Elongation_max = LinkPropertiesStruct.MaxLinkElongation;

% number of decoupled points
Ndecoupled = Noutput*2;

% initialize blank force array
BlankForceArray = [(LatticeGeometryStruct.inputNodes)' zeros(length(LatticeGeometryStruct.inputNodes),6)];

% compute MSE for each decoupled point
for decoupledIDX = 1:Ndecoupled

    % set single force for MSE computation, all zeroes but one
    if mod(decoupledIDX,2) == 0
        DecoupledBehaviorStruct.forces = BlankForceArray;
        DecoupledBehaviorStruct.forces((decoupledIDX/2),3) = 1;
    else
        DecoupledBehaviorStruct.forces = BlankForceArray;
        DecoupledBehaviorStruct.forces(((decoupledIDX+1)/2),2) = 1;
    end

    % dummy target 
    DecoupledBehaviorStruct.Target = DenseSweepStruct.Target(:,:,1);

    % package FEM
    [DOFFinal,Final, F, DOFnodes,k_base,R6,RT6,DOF,Degrees_per_element] = FEM_setup(LinkPropertiesStruct,LatticeGeometryStruct,DecoupledBehaviorStruct);
    FEMStructVia.DOFFinal=DOFFinal;
    FEMStructVia.Final=Final;
    FEMStructVia.F=F;
    FEMStructVia.DOFnodes=DOFnodes;
    FEMStructVia.k_base=k_base;
    FEMStructVia.R6=R6;
    FEMStructVia.RT6=RT6;
    FEMStructVia.DOF=DOF;
    FEMStructVia.Degrees_per_element=Degrees_per_element;

     % compute node displacements for discrete point
    deformedDecoupled(:,:,:,decoupledIDX)=FEM(LinkPropertiesStruct, LatticeGeometryStruct, DecoupledBehaviorStruct,FEMStructVia,OptimizerDataStruct,xFinal);

end

% extract coordinates of deformed output nodes
OutDeformedDecoupled = deformedDecoupled(LatticeGeometryStruct.outputNodes,:,1,:);

% compute displacement of each output node
for outIDX = 1:Noutput
    % compute displacements
    for decoupledIDX = 1:Ndecoupled
        OutDispDecoupled(outIDX,:,1,decoupledIDX) = OutDeformedDecoupled(outIDX,:,1,decoupledIDX)-OutCoordInitial(outIDX,:);
    end
end

% generate reduced stiffness matrix
% U = K * x
% U column vector [Ux1 Uy1 Ux2 Uy2 ...] where 1, 2, etc is node number
% F column vector [Fx1 Fy1 Fx2 Fy2 ...] where 1, 2, etc is node number
for decoupledIDX = 1:Ndecoupled
    for outIDX = 1:Noutput
        Kf2u((outIDX*2-1):(outIDX*2),decoupledIDX) = [OutDispDecoupled( outIDX,1,1,decoupledIDX ); OutDispDecoupled( outIDX,2,1,decoupledIDX )];
    end
end

%% Validation & Singular Value Decomposition

% singular value decomposition
[U,S,V] = svd(Kf2u);

% number of points in circle
Ncirc = 1000;

% set of input angles
theta = linspace(0,2*pi,Ncirc);

% compute individual transformations and final
% A = U*S*V'
for circIDX = 1:Ncirc
    % compute force magnitudes
    Ymag = cos( theta(circIDX) );
    Xmag = sin( theta(circIDX) );

    % arrange force input vector
    F = [];
    for outIDX = 1:Noutput
        F = [F;Xmag;Ymag];
    end
    Farray(:,circIDX) = F;

    % compute output displacement
    u = Kf2u*F;
    Uarray(:,circIDX) = u;
    
    % compute intermediate transformations from SVD
    trans(:,circIDX,1) = F;
    trans(:,circIDX,2) = (V')*F;
    trans(:,circIDX,3) = S*(V')*F;
    trans(:,circIDX,4) = U*S*(V')*F;
end

figure
hold on
for outIDX = 1:Noutput
    plot(Uarray( (outIDX*2-1) ,:), Uarray((outIDX*2) ,:),'b-')
end

figure
hold on
for outIDX = 1:Noutput
    for transIDX = 1:4
        subplot(Noutput,4, (outIDX-1)*4+(transIDX) )
        plot(trans( (outIDX*2-1) ,:,transIDX), trans((outIDX*2) ,:,transIDX),'r-','LineWidth',1.5)
        axis equal
        box on
        set(gca,'LineWidth',1.5)
        set(gca,'XTick',[])
        set(gca,'YTick',[])
    end
end

%% Plot Elliptical Space

% input force magnitude space
FmagArray = linspace(0,2,5);

% compute Force -> Displkacement transformation for a full revolution at
% each magnitude
for FIDX = 1:length(FmagArray)
    for circIDX = 1:Ncirc

        % compute force magnitudes
        Ymag = cos( theta(circIDX) ) * FmagArray(FIDX);
        Xmag = sin( theta(circIDX) ) * FmagArray(FIDX);

        % arrange force input vector
        F = [];
        for outIDX = 1:Noutput
            F = [F;Xmag;Ymag];
        end
        Flinspace(:,circIDX,FIDX) = F;

        % compute output displacement
        u = Kf2u*F;
        Ulinspace(:,circIDX,FIDX) = u;
    end
end

% color scale for space
Color0 = [0 92 105] ./ 255;
Colorf = [92 0 105] ./ 255;


% colored plot
figure
hold on
for FIDX = 1:length(FmagArray)
    % select color 
    colorg = Color0 + ((Colorf - Color0) / (max(FmagArray)-min(FmagArray)) * (FmagArray(FIDX)-min(FmagArray)));

    for outIDX = 1:Noutput
        % extract displacements
        dispX = Ulinspace( (outIDX*2-1) ,:,FIDX);
        dispY = Ulinspace( (outIDX*2) ,:,FIDX);

        % add displacements to original positions
        posX = OutCoordInitial(outIDX,1) + dispX;
        posY = OutCoordInitial(outIDX,2) + dispY;

        % add to subplot
        subplot(Noutput,1, Noutput+1-outIDX )
        hold on
        plot(posX,posY,"Color",colorg)
        grid on
        axis([mean(BehaviorStruct.Target(outIDX,1,:))-0.002 mean(BehaviorStruct.Target(outIDX,1,:))+0.002 BehaviorStruct.Target(outIDX,2,1)-0.0005 BehaviorStruct.Target(outIDX,2,1)+0.0005])
    end
end

%% plot singular values

figure
hold on

for singIDX = 1:size(S,1)
    plot(singIDX,S(singIDX,singIDX),'*k','MarkerSize',15,'LineWidth',2.5)
end
box on
set(gca,'LineWidth',2.5)
yscale('log')
set(gca,'FontSize',18)


%% Evaluate Reduced Matrix Properties

% distribution of input effects
WeightingVector = U\diag(S);
figure
plot(abs(WeightingVector)/sum(abs(WeightingVector)),'*k')


