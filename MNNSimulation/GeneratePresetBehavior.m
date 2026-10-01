clear
clc
close all

% Author: Pietro Sainaghi
% code modified from: Binary_MinStiffness by Ryan Lee
% works with the other functions in this folder, written by Reinier Kuppens
% and Ryan Lee

% this code studies a given MNN lattice ability to learn behaviors with
% respect to various lattice parameters

%% Add Function Folders

addpath(genpath([pwd,'\Main_Functions']))

%% Configure settings for created behaviors

% number of sets behaviors to create
Nfiles = 1;

% type of RNG for randomized behavior
RNGtype = 'deterministic'; % options: 'dateandtime' 'deterministic' 'nocontrol'

% assemble PregenBehStruct
PregenBehStruct.Nfiles = Nfiles;
PregenBehStruct.RNGtype = RNGtype;



%% Individual Link Parameters

% these should be set based on CADs and FEA for the links tobe used in the
% MNN

% length of each beam
L1                  = 0.156; % m

% assemble LinkPropertiesStruct
LinkPropertiesStruct.L1 = L1;


% assemble PregenBehStruct
PregenBehStruct.L1 = L1;



%% Lattice Gemetry Parameters

% select type of lattice to simulate 
    % only implemented for triangular lattice at the moment
latticeType          = 1; %1) triangular 2) square

% number of nodes that receive forces and measure displacements
    % set as array to test multiple ones at once
NinputANDoutput     = 8;

% array of numbers of layers
    % set as array to test multiple ones at once
Nlayers              = 8;

% dimensions in space
DOI = 3;

% assemble LatticeGeometryStruct
LatticeGeometryStruct.latticeType = latticeType;
LatticeGeometryStruct.NinputANDoutput = NinputANDoutput;
LatticeGeometryStruct.Nlayers = Nlayers;
LatticeGeometryStruct.DOI = DOI;

% assemble PregenBehStruct
PregenBehStruct.latticeType = latticeType;
PregenBehStruct.NinputANDoutput = NinputANDoutput;
PregenBehStruct.Nlayers = Nlayers;
PregenBehStruct.DOI = DOI;


%% Behavior Configuration

% type of study 
caseType        = 4;  % 1) sinusoid behavior     2) random forces   4) generalized

% number of force behaviors in randomized behaviors 
% number of discrete via points in generalized curve
Ncases             = 1100; 

% minimum force allowed as fraction of maximum force
threshold       = 0.3;

% maximum allowed elongation of random behavior
    % set as array to test multiple ones at once
Elongation_max  = 0.0025; % units of m

% aplitude of sine wave for sinusoid behavior
    % set same value as Elongation_max for looping convenience, since
    % random and sinusoid are not tested in the same run
dx              = Elongation_max; % units of m

% maximum allowed input force
    % set as array to test multiple ones at once
MaxForce = 1*8; % units of N

% assemble BehaviorStruct
BehaviorStruct.caseType = caseType;
BehaviorStruct.Ncases = Ncases;
BehaviorStruct.threshold = threshold;
BehaviorStruct.Elongation_max = Elongation_max;
BehaviorStruct.dx = dx;
BehaviorStruct.MaxForce = MaxForce;
BehaviorStruct.RNGtype = RNGtype;

% assemble PregenBehStruct
PregenBehStruct.caseType = caseType;
PregenBehStruct.Ncases = Ncases;
PregenBehStruct.threshold = threshold;
PregenBehStruct.Elongation_max = Elongation_max;
PregenBehStruct.dx = dx;
PregenBehStruct.MaxForce = MaxForce;






%% select folder to save .mat files

if PregenBehStruct.caseType == 1
    fileName = 'SinusoidBehavior';
elseif PregenBehStruct.caseType == 2
    fileName = 'randomBehavior';
    switch PregenBehStruct.RNGtype
        case 'dateandtime'
            fileName = [fileName,'_',num2str(Ncases),'beh_datetimerng'];
        case 'deterministic'
            fileName = [fileName,'_',num2str(Ncases),'beh_deterministicrng'];
        case 'nocontrol'
            fileName = [fileName,'_',num2str(Ncases),'beh_nocontrol'];
    end
elseif PregenBehStruct.caseType == 4
    fileName = 'generalizedBehavior';
    fileName = [fileName,'_',num2str(Ncases),'viapoints'];
end
fileType = '.mat';
t= datetime;
t.Format = 'yyyy-MM-dd_@_HHmm-ss';
t = char(t);
loopName = [fileName,'_','run@',t] ;

savePath = uigetdir();


%% generate behavior files

for fileIDX = 1:Nfiles
    
    % set deterministic RNG for random behavior
    switch BehaviorStruct.RNGtype
        case 'deterministic'
            BehaviorStruct.RNG_seed_det = fileIDX*3;
    end
    % generate MNN lattice structure
    LatticeGeometryStruct = packageLattice(LinkPropertiesStruct, LatticeGeometryStruct);
    % number of nodes in lattice
    PregenBehStruct.Ncoord = LatticeGeometryStruct.Ncoord;
    % number of beams in the lattice
    PregenBehStruct.Nbeams = LatticeGeometryStruct.Nbeams;
    % coordinates of each node when lattice is undisturbed
    PregenBehStruct.coord_initial = LatticeGeometryStruct.coord_initial;
    % list of nodes each link connects
    PregenBehStruct.connectivity = LatticeGeometryStruct.connectivity;
    % indices of ground nodes
    PregenBehStruct.bound = LatticeGeometryStruct.bound;
    % indices of input nodes
    PregenBehStruct.inputNodes = LatticeGeometryStruct.inputNodes;
    % indices of output nodes
    PregenBehStruct.outputNodes = LatticeGeometryStruct.outputNodes;
    
    % generate behavior
    BehaviorStruct = packageBehavior(LatticeGeometryStruct,BehaviorStruct);
    % target node coordinates
    PregenBehStruct.Target = BehaviorStruct.Target;
    % forces locations in DOF vector format
    PregenBehStruct.forces = BehaviorStruct.forces;
    
    % format file name
    fullName = [savePath,'\',loopName,'_',num2str(fileIDX),fileType]
    
    % save to file
    save(fullName,'PregenBehStruct')
    
end

