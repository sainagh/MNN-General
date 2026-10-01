clear
clc
close all

% Author: Pietro Sainaghi
% code modified from: Binary_MinStiffness by Ryan Lee
% works with the other functions in this folder, mostly written by Ryan Lee

% this code studies a given MNN lattice ability to learn behaviors with
% respect to various lattice parameters

%% Add Function Folders

addpath(genpath([pwd,'\Main_Functions']))

%% Individual Link Parameters

% these should be set based on CADs and FEA for the links tobe used in the
% MNN

% length of each beam
L1                  = 0.011; % m
%X-component of beam under 60 degrees / 1/3 pi rad
L2                  = L1*cosd(60); 
%Y-component of beam under 60 degrees / 1/3 pi rad
L3                  = L1*cosd(30);     

% maximum link elongation
MaxLinkElongation   = 10e-6; % m

% maximum link compression
MaxLinkCompression   = 10e-6; % m

% stiffness values
    % Force acting perpendicular while rotation is fixed: utcross=0;
k55                  = 24.57366; %N/m
    % Force acting perpendicular with rotation of tip allowed
k56                  = 0.120876; %N/rad
    % Moment acting on end while perpendicular translation of tip is allowed 
k66                  = 41.57994; %N*m/rad

% passive axial stiffness
kLinPassive             = 2781.602; %N/m

% maximum controllable stiffness
kLinMax                 = 3000; %N/m

% minimum controllable stiffness
kLinMin                 = -2000; %N/m

% assemble LinkPropertiesStruct
LinkPropertiesStruct.L1 = L1;
LinkPropertiesStruct.L2 = L2;
LinkPropertiesStruct.L3 = L3;
LinkPropertiesStruct.MaxLinkElongation = MaxLinkElongation;
LinkPropertiesStruct.MaxLinkCompression = MaxLinkCompression;
LinkPropertiesStruct.k55 = k55;
LinkPropertiesStruct.k56 = k56;
LinkPropertiesStruct.k66 = k66;
LinkPropertiesStruct.kLinPassive = kLinPassive;
LinkPropertiesStruct.kLinMax = kLinMax;
LinkPropertiesStruct.kLinMin = kLinMin;

%% Nonlinearities in Beam Axial Stiffness

% wether to have nonlinear stiffness
nonlinearStiffness = false;

% type of nonlinear function for beam stiffness
    % 'ATAN' Arctangent: F = Klin * u + Knonlin * tan(u)
    % 'RAMP' Ramp function: F = Klin * u + Knonlin * u * (u>=0)
    % 'EXP' Exponential: F = Klin * u + Knonlin * Fmax * (exp(u)-1)
    % 'CUBE' Cubic: F = Klin * u + Knonlin * Fmax * u^3
    % WIP 'ERF' Error function: F = K * ( 2/sqrt(pi) * integral( e ^ -( x / MaxLinkElongation ) ^ 2 ) )
    % WIP 'CB' Cubic: uses F = K * u^3;
nonLinearityType = 'ATAN';

% assemble LinkPropertiesStruct
LinkPropertiesStruct.nonlinearStiffness = nonlinearStiffness;
LinkPropertiesStruct.nonLinearityType = nonLinearityType;

%% Lattice Gemetry Parameters

% select type of lattice to simulate 
latticeType          = 1; %1) triangular 2) square

% number of nodes that receive forces and measure displacements
    % set as array to test multiple ones at once
NinputANDoutputArray      = [2];

% array of numbers of layers
    % set as array to test multiple ones at once
NlayersArray              = [2];

% dimensions in space
DOI = 3;

% assemble LatticeGeometryStruct
LatticeGeometryStruct.latticeType = latticeType;
LatticeGeometryStruct.NinputANDoutputArray = NinputANDoutputArray;
LatticeGeometryStruct.NlayersArray = NlayersArray;
LatticeGeometryStruct.DOI = DOI;

%% Behavior Configuration

% type of study 
caseType        = 1;  % 1) sinusoid behavior     2) random forces       3) saved behavior
% GeneratePresetBehavior will create .mat files compatible with caseType 3

% number of force behaviors in randomized behaviors 
    % set as array to test multiple ones at once
    % can only be 2 for caseType = 1
    % will be overwritten for caseType = 3
NcasesArray             = [2]; 

% minimum force allowed as fraction of maximum force
threshold       = 0.3;

% maximum allowed elongation of random behavior
    % set as array to test multiple ones at once
    % will be overwritten for caseType = 3
Elongation_maxArray  = [2.5e-6]; % units of m

% aplitude of sine wave for sinusoid behavior
    % set same value as Elongation_maxArray for looping convenience, since
    % random and sinusoid are not tested in the same run
    % will be overwritten for caseType = 3
dxArray              = Elongation_maxArray; % units of m

% maximum allowed input force
    % set as array to test multiple ones at once
    % will be overwritten for caseType = 3
MaxForceArray = [1e-3]; % units of N

% type of RNG for randomized behavior
RNGtype = 'dateandtime'; % options: 'dateandtime' 'deterministic' 'nocontrol'


% assemble BehaviorStruct
BehaviorStruct.caseType = caseType;
BehaviorStruct.NcasesArray = NcasesArray;
BehaviorStruct.threshold = threshold;
BehaviorStruct.Elongation_maxArray = Elongation_maxArray;
BehaviorStruct.dxArray = dxArray;
BehaviorStruct.MaxForceArray = MaxForceArray;
BehaviorStruct.RNGtype = RNGtype;



%% Optimizer Configuration

% force scaling
ForceScaling        = false;

% make a physically possible network
    % max elongation in each beam will be contrained to be less than Elongation_maxArray
    % some optimizers don't work well with the constraint setup
EnforceMaxElongation = false;

% number of runs for each set of behaviors
startPts            = 1;

% type of initial condition
icType          = 3; %1) all max, 2) all min, 3) rand

% whether to use deterministic random initial conditions
    % use date and time if you want to ensure different initial conditions every time
    % use seed multiplier if you want to control your initial condition     
deterministicRNG = 2; % 0) use base RNG 1) specify a seed for rng 2) seed is date and time
RNGseedMult = 3; % for deterministicRNG = 1, the seed will be loopnumber*RNGseedMult

% wether to use MSE or scaled normalized MSE
doScaleMSE = 0; % 0) use MSE, 1) use unitless MSE

% select units for MSE
    % 1 = m, 100 = cm, 1000 = mm, 1000000 = um
MSEunits = 1000000;

% whether to run a discrete or continuous optimization
    % discrete only works with GA
discORcont = 'discrete'; % 'discrete' or 'continuous'

% set diiscrete set of stiffness values
    % unused if discORcont = 'continuous'
    % set as cell to sweep for different stiffness libraries
PossibleStiffnessCell{1} = linspace(kLinMin, kLinMax, 5000); %N/m


% type of optimizer to use
    % set as cell to test multiple ones at once
    % all performance evaluations depend on hyper-parameter tuning
    % GA: genetic algorithm, slowest but very accurate
    % SQP: sequential quadratic programming, fast but less inaccurate
    % FPS: pattern search using matlab function, fast and medium accuracy
    % PPS: partial pattern search, developed by Ryan H. Lee, slow but accurate TODO NYI
    % AGD: analytical gradient method, developed by Jiaji Chen, fast
optimizerArray = {'GA'}; % 'GA', 'SQP', 'FPS', 'PPS', 'AGD'

% Genetic Algorithm Hyper-Parameters
% Only used if optimizer is GA
GApopulation = 250; % number of function evaluations each iteration
GAgenerations = 10000; % max number of iterations
GAstallgenerations = 100; % max number of iterations with same outcome for exit flag
GAparallelPool = false; % whether to use multithreaded optimizer (faster but more resource intensive)
GAerrorChangeThreshold  = 1e-100; % change in error to indicate convergence

% Sequential Quadratic Programming Hyper-Parameters
% Only used if optimizer is SQP
% TODO add other optimization parameters as well
SQPmaxIterations = 2000000; % number of optimizer iterations
SQPmaxFunEvals = 1e8; % maximum number of allowed function evaluations
SQPerrorChangeThreshold = 1e-16; % change in error to indicate convergence
SQPparallelization = true; % whether to use multithreaded optimizer (faster but more resource intensive)
SQPsteptol = 1e-16; % change in step size to indicate convergence

% Full Pattern Search Hyper-Parameters
FPSerrorChangeThreshold = 1e-6; % change in error to indicate convergence
FPSmaxIterations = 1e6; % maximum iterations
FPSsearchFunc = []; % possible pattern types for a given iteration, options: "GPSPositiveBasis2N" | "GPSPositiveBasisNp1" | "GSSPositiveBasis2N" | "GSSPositiveBasisNp1" | "MADSPositiveBasis2N" | "MADSPositiveBasisNp1" | "searchga" | "searchlhs" | "searchneldermead" | "rbfsurrogate" | {[]} |
FPSalgorithm = 'classic'; % pattern search algorithms, options: {"classic"} | "nups" | "nups-gps" | "nups-mads"
FPSParallelPool = false; % whether to use multithreaded optimizer (faster but more resource intensive)

% Analytical Gradient Method Hyper-Parameters
AGDmaxIterations = 1000;

% whether to set computational cache in separate folder
% will be same as results if false
% set to true if using GA or other algorithm with lots of function
% evaluations
SpecifyCache = 0;


% assemble OptimizerDataStruct
OptimizerDataStruct.ForceScaling = ForceScaling;
OptimizerDataStruct.EnforceMaxElongation = EnforceMaxElongation;
OptimizerDataStruct.startPts = startPts;
OptimizerDataStruct.icType = icType;
OptimizerDataStruct.doScaleMSE = doScaleMSE;
OptimizerDataStruct.discORcont = discORcont;
OptimizerDataStruct.optimizerArray = optimizerArray;
OptimizerDataStruct.deterministicRNG = deterministicRNG;
OptimizerDataStruct.PossibleStiffnessCell = PossibleStiffnessCell;
OptimizerDataStruct.MSEunits = MSEunits;
OptimizerDataStruct.SpecifyCache = SpecifyCache;
OptimizerDataStruct.GAerrorChangeThreshold = GAerrorChangeThreshold;
OptimizerDataStruct.GApopulation = GApopulation;
OptimizerDataStruct.GAgenerations = GAgenerations;
OptimizerDataStruct.GAstallgenerations = GAstallgenerations;
OptimizerDataStruct.GAparallelPool = GAparallelPool;
OptimizerDataStruct.RNGseedMult = RNGseedMult;
OptimizerDataStruct.SQPmaxIterations = SQPmaxIterations;
OptimizerDataStruct.SQPmaxFunEvals = SQPmaxFunEvals;
OptimizerDataStruct.SQPerrorChangeThreshold = SQPerrorChangeThreshold;
OptimizerDataStruct.FPSerrorChangeThreshold = FPSerrorChangeThreshold;
OptimizerDataStruct.FPSmaxIterations = FPSmaxIterations;
OptimizerDataStruct.FPSsearchFunc = FPSsearchFunc;
OptimizerDataStruct.FPSalgorithm = FPSalgorithm;
OptimizerDataStruct.FPSParallelPool = FPSParallelPool;
OptimizerDataStruct.SQPparallelization = SQPparallelization;
OptimizerDataStruct.SQPsteptol = SQPsteptol;

%% Plotting Options Configuration

% scale width of deformed parts
    % Not yet implemented
scalingterm = 0;

% plot undeformed lattice
plotUndeformed = 1;

% plot stiffness combinations in lattice
plotColoredLattice = 1;

% plot deformed lattice
    % set to 1 only if the loop inludes one lattice and one set of behaviors
plotDeformed = 0;

% plot endpoints
    % set to 1 only if the loop inludes lattice and one set of behaviors
plotEndpoints = 1;
    % amplitude of axes around initial position
plotEndPointsAmplitude = 5e-6; % in m

% assemble plotOptionsStruct
plotOptionsStruct.scalingterm = scalingterm;
plotOptionsStruct.plotUndeformed = plotUndeformed;
plotOptionsStruct.plotDeformed = plotDeformed;
plotOptionsStruct.plotEndpoints = plotEndpoints;
plotOptionsStruct.plotColoredLattice = plotColoredLattice;
plotOptionsStruct.plotEndPointsAmplitude = plotEndPointsAmplitude;


%% Plot Nonlinearity Within Viable Range


if nonlinearStiffness
    figure(999)
    hold on
    plot([-max(MaxLinkElongation,MaxLinkCompression) max(MaxLinkElongation,MaxLinkCompression)], [-max(MaxForceArray) max(MaxForceArray)],'b--')
    plot([-max(MaxLinkElongation,MaxLinkCompression) -max(MaxLinkElongation,MaxLinkCompression)], [-max(MaxForceArray) max(MaxForceArray)],'k-')
    plot([max(MaxLinkElongation,MaxLinkCompression) max(MaxLinkElongation,MaxLinkCompression)], [-max(MaxForceArray) max(MaxForceArray)],'k-')
    plot([-max(MaxLinkElongation,MaxLinkCompression) max(MaxLinkElongation,MaxLinkCompression)], [-max(MaxForceArray) -max(MaxForceArray)],'k-')
    plot([-max(MaxLinkElongation,MaxLinkCompression) max(MaxLinkElongation,MaxLinkCompression)], [max(MaxForceArray) max(MaxForceArray)],'k-')
    u_plot = linspace(-max(MaxLinkElongation,MaxLinkCompression)*1.1,max(MaxLinkElongation,MaxLinkCompression)*1.1,500);
    switch nonLinearityType
        case 'ATAN'
            F_plot = 2 * max(MaxForceArray) / pi * atan(u_plot/max(MaxLinkElongation,MaxLinkCompression)*pi/2);
        case 'EXP'
            F_plot = max(MaxForceArray) * (exp(u_plot/max(MaxLinkElongation,MaxLinkCompression))-1);
        case 'CUBE'
            F_plot = max(MaxForceArray) * (u_plot/max(MaxLinkElongation,MaxLinkCompression)).^3;
        case 'RAMP'
            F_plot = [zeros(1,250) linspace(0,max(MaxForceArray)*1.1,250)];
    end
    plot(u_plot, F_plot,'r-')
    grid on
    axis([-max(MaxLinkElongation,MaxLinkCompression)*1.1 max(MaxLinkElongation,MaxLinkCompression)*1.1 -max(MaxForceArray)*1.1 max(MaxForceArray)*1.1])
end


%% Save File

% flag to save learning output
IWantToSaveOutput = true;

% flag to autoselect output
    AutoSpecifyOutput = false;

if IWantToSaveOutput
    
    % ouput location is specified automatically
    if AutoSpecifyOutput
        fileName = 'PlanarMeso_learningAbilityStudy';
        fileType = '.mat';
        t= datetime;
        t.Format = 'yyyy-MM-dd_@_HHmm-ss';
        t = char(t);
        fullName = [fileName,'_','run@',t, fileType] ;

        savePath = [pwd,'\Results'];
        folderName = 'MNN_study';
        fullPath = [savePath,'\',folderName,'\', fullName];

        mkdir([savePath,'\',folderName]);

    % output location is specified manually
    else
        fileName = 'PlanarMeso_learningAbilityStudy';
        fileType = '.mat';
        t= datetime;
        t.Format = 'yyyy-MM-dd_@_HHmm-ss';
        t = char(t);
        fullName = [fileName,'_','run@',t, fileType] ;

        savePath = uigetdir(pwd,'Select output save location');
        folderName = 'MNN_study';
        fullPath = [savePath,'\',folderName,'\', fullName];

        mkdir([savePath,'\',folderName]);
    end

end

%% select pregen behaviors

AutoSelectBeh = false;
AutoSelectedBehString = [pwd,'\SavedBehaviors\generalizedBehavior_3viapoints_run@2023-09-22_@_1448-31_1.mat'];

% only happens if caseType = 3
if caseType == 3

    if AutoSelectBeh
        load(AutoSelectedBehString);
        % BehaviorStruct.fileList = fileList;
        % BehaviorStruct.dataPath = dataPath;
        NcasesArray = [PregenBehStruct.Ncases];
        BehaviorStruct.NcasesArray = [PregenBehStruct.Ncases];
        NumBehFiles = 1;
    else
        % TODO
        % To improve
        [fileList, dataPath] = uigetfile('*.mat','Select preset behavior','MultiSelect', 'on');
        NumBehFiles = 1;
        load([dataPath,fileList]);
        BehaviorStruct.fileList = fileList;
        BehaviorStruct.dataPath = dataPath;
        NcasesArray = [PregenBehStruct.Ncases];
        BehaviorStruct.NcasesArray = [PregenBehStruct.Ncases];
    end
    
    % check if preset behavior and lattice geometry setting set in code
    % here are compatible
    if (length(LatticeGeometryStruct.NinputANDoutputArray) > 1) || (length(LatticeGeometryStruct.NlayersArray) > 1)
        error('You must only specify one lattice size if you are using preset behaviors')
    end
    if LatticeGeometryStruct.NinputANDoutputArray ~= PregenBehStruct.NinputANDoutput
        error('Preset behavior structure is incompatible with your current lattice geometry, wrong NinputANDoutput')
    end
    if LatticeGeometryStruct.NlayersArray ~= PregenBehStruct.Nlayers
        error('Preset behavior structure is incompatible with your current lattice geometry, wrong Nlayers')
    end
else
    NumBehFiles = 1;
end

%% Select Computation Cache

if OptimizerDataStruct.SpecifyCache == 0
    
    fileNameCache = 'PlanarMeso_learningAbilityStudy';
    fileType = '.mat';
    t= datetime;
    t.Format = 'yyyy-MM-dd_@_HHmm-ss';
    t = char(t);
    cacheName = [fileNameCache,'_','run@',t, '_chachedata',fileType] ;
    
    savePathCache = savePath;
    folderNameCache = 'MNN_study_cache';
    
    CachePath = [savePathCache,'\',folderNameCache,'\', cacheName];
    
    mkdir([savePathCache,'\',folderNameCache]);
    
elseif OptimizerDataStruct.SpecifyCache == 1
    
    fileNameCache = 'PlanarMeso_learningAbilityStudy';
    fileType = '.mat';
    t= datetime;
    t.Format = 'yyyy-MM-dd_@_HHmm-ss';
    t = char(t);
    cacheName = [fileNameCache,'_','run@',t, '_chachedata',fileType] ;
    
    savePathCache = uigetdir(pwd,'Select cache location');
    folderNameCache = 'MNN_study_cache';
    
    CachePath = [savePathCache,'\',folderNameCache,'\', cacheName];
    
    mkdir([savePathCache,'\',folderNameCache]);
end



save('CacheDirectory.mat','CachePath');


%% Main Loops

% overall run counter
NoptRuns = 0;

% loop through number of inputs and output nodes
for NinputANDoutputIter = 1:length(NinputANDoutputArray)
    LatticeGeometryStruct.NinputANDoutput = LatticeGeometryStruct.NinputANDoutputArray(NinputANDoutputIter);

    % loop through number of layers
    for NlayersIter = 1:length(NlayersArray)
        LatticeGeometryStruct.Nlayers = LatticeGeometryStruct.NlayersArray(NlayersIter);

        % loop through various stiffness libraries
        for StiffnessLibraryIter = 1:length(PossibleStiffnessCell)
            OptimizerDataStruct.PossibleStiffnessArray = OptimizerDataStruct.PossibleStiffnessCell{StiffnessLibraryIter};

            % loop through output amplitudes
            for dxIter = 1:length(dxArray)
                BehaviorStruct.dx = BehaviorStruct.dxArray(dxIter);
                BehaviorStruct.Elongation_max = BehaviorStruct.Elongation_maxArray(dxIter);

                % loop through input force magnitudes
                for MaxForceIter = 1:length(MaxForceArray)
                    BehaviorStruct.MaxForce = BehaviorStruct.MaxForceArray(MaxForceIter);

                    % loop through number of cases
                    for NcasesIter = 1:length(NcasesArray)
                        BehaviorStruct.Ncases = BehaviorStruct.NcasesArray(NcasesIter);

                        % loop though optimizers
                        for optimizerIter = 1:length(OptimizerDataStruct.optimizerArray)
                            OptimizerDataStruct.optimizer = OptimizerDataStruct.optimizerArray{optimizerIter};

                            % Loop through number of behaviors
                            for behaviorSetIter = 1:NumBehFiles

                                % load behavior if from file
                                if caseType == 3
                                    if AutoSelectBeh
                                        load(AutoSelectedBehString);
                                    else
                                        load([dataPath,fileList])
                                    end
                                    BehaviorStruct.Target = PregenBehStruct.Target;
                                    BehaviorStruct.forces = PregenBehStruct.forces;
                                    BehaviorStruct.Ncases = PregenBehStruct.Ncases;
                                    BehaviorStruct.Elongation_max = PregenBehStruct.Elongation_max;
                                    BehaviorStruct.MaxForce = PregenBehStruct.MaxForce;
                                    BehaviorStruct.dx = PregenBehStruct.dx;
                                end

                                % set overall run number
                                NoptRuns = NoptRuns + 1;
                                % debug study display
                                disp(['Design study run number: ', num2str(NoptRuns),' / ',num2str(length(optimizerArray)*length(NcasesArray)*length(MaxForceArray)*length(dxArray)*length(NlayersArray)*length(NinputANDoutputArray))])
                                disp(['Number of input/output nodes: ', num2str(LatticeGeometryStruct.NinputANDoutput)])
                                disp(['Number of layers: ', num2str(LatticeGeometryStruct.Nlayers)])
                                disp(['Stiffness Library Number: ', num2str(StiffnessLibraryIter), ' / ', num2str(length(OptimizerDataStruct.PossibleStiffnessCell))])
                                if BehaviorStruct.caseType == 1
                                    disp(['Behavior type: sinusoid'])
                                elseif BehaviorStruct.caseType == 2
                                    disp(['Behavior type: random',', ',BehaviorStruct.RNGtype, ' RNG'])
                                elseif BehaviorStruct.caseType == 3
                                    disp(['Behavior type: saved'])
                                end
                                disp(['Number of sets of behaviors: ', num2str(NumBehFiles)])
                                disp(['Behavior Output Amplitude: ', num2str(BehaviorStruct.Elongation_max)])
                                disp(['Behavior max force: ', num2str(BehaviorStruct.MaxForce)])
                                disp(['Number of behaviors to optimize: ', num2str(BehaviorStruct.Ncases)])
                                disp(['Type of optimizer: ', (OptimizerDataStruct.optimizer)])
                                disp(['Optimization space: ', (OptimizerDataStruct.discORcont)])

                                % set deterministic RNG for random behavior
                                switch BehaviorStruct.RNGtype
                                    case 'deterministic'
                                        BehaviorStruct.RNG_seed_det = NoptRuns*3;
                                end


                                % generate MNN lattice structure
                                LatticeGeometryStruct = packageLattice(LinkPropertiesStruct, LatticeGeometryStruct);


                                if (caseType == 1) || (caseType == 2)
                                    % generate behavior
                                    BehaviorStruct = packageBehavior(LatticeGeometryStruct,BehaviorStruct);
                                end



                                % compute MSE scale factor
                                % unused if OptimizerDataStruct.doScaleMSE = 0
                                if OptimizerDataStruct.doScaleMSE == 1
                                    MSEscale = computeMSEScaleFactor(LatticeGeometryStruct,BehaviorStruct);
                                    % store in structure OptimizerDataStruct
                                    OptimizerDataStruct.MSEscale = MSEscale;
                                else
                                    % store in structure OptimizerDataStruct
                                    % zero means unused
                                    OptimizerDataStruct.MSEscale = 0;
                                end



                                % set up finite truss structure (finite element used equivalently)
                                [DOFFinal,Final, F, DOFnodes,k_base,R6,RT6,DOF,Degrees_per_element] = FEM_setup(LinkPropertiesStruct,LatticeGeometryStruct,BehaviorStruct);

                                %store in structure FEMStruct
                                FEMStruct.DOFFinal = DOFFinal;
                                FEMStruct.Final = Final;
                                FEMStruct.F = F;
                                FEMStruct.DOFnodes = DOFnodes;
                                FEMStruct.k_base = k_base;
                                FEMStruct.R6 = R6;
                                FEMStruct.RT6 = RT6;
                                FEMStruct.DOF = DOF;
                                FEMStruct.Degrees_per_element = Degrees_per_element;

                                % run MNN with given specifications
                                ResultsStruct = MNN_run(...
                                    LinkPropertiesStruct,...
                                    LatticeGeometryStruct,...
                                    BehaviorStruct,...
                                    OptimizerDataStruct,...
                                    FEMStruct,...
                                    plotOptionsStruct);

                                % save to individual cell
                                SweepOutput{...
                                    NinputANDoutputIter,...
                                    NlayersIter,...
                                    StiffnessLibraryIter,...
                                    dxIter,...
                                    MaxForceIter,...
                                    NcasesIter,...
                                    optimizerIter} = ResultsStruct;


                                disp(newline)
                                disp(newline)

                            end % end behavior set loop
                        end % end optimizer loop
                    end % end number of behaviors loop
                end % end input force loop
            end % end output amplitude loop
        end % end stiffness libraries loop
    end % end number of layers loop
end % end number of input nodes loop

%% Save Ouput to File

if IWantToSaveOutput
    
    save(fullPath,'LinkPropertiesStruct','LatticeGeometryStruct','BehaviorStruct','OptimizerDataStruct','FEMStruct','plotOptionsStruct','SweepOutput')
end














