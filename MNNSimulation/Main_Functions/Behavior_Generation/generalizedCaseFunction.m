function [Target, forces] = generalizedCaseFunction...
        (LatticeGeometryStruct,BehaviorStruct)
    
% function written by Pietro Sainaghi
% generates behavior set that constritutes a set of via points in a
% generalized force-amplitude space
    
%% function inputs

% lattice structure - LatticeGeometryStruct
inputNodes = LatticeGeometryStruct.inputNodes;
outputNodes = LatticeGeometryStruct.outputNodes;
nLayers = LatticeGeometryStruct.Nlayers;
DOI = LatticeGeometryStruct.DOI;
coord_initial = LatticeGeometryStruct.coord_initial;
nCoord = LatticeGeometryStruct.Ncoord;

% behavior structure - BehaviorStruct
nCases = BehaviorStruct.Ncases;
dx = BehaviorStruct.dx;
maxForce = BehaviorStruct.MaxForce;



%% formatting

nIO = length(inputNodes); % number of inputs and outputs
forces = zeros(nIO,7,nCases);

%% generalized curve

% start and end of generalized curve
startGEN = dx;
finishGEN = -dx;

% magnitude of rotating force
Fmag=maxForce/nIO;

% theta via points
thetaVIA = linspace(0,pi,nCases);

% amplitude via points
ampVIA = dx * cos( thetaVIA );

figure
plot(thetaVIA, ampVIA,'d-')
axis([-0.2 pi+0.2 -1.5*dx 1.5*dx])
xlabel('Force Angle [rad]')
ylabel('Output Displacement Amplitude [m]')
grid on


%% set input forces

% populate array for each via point
for caseIDX = 1:nCases
    % set indices of nodes receiving forces
    forces(:,1,caseIDX) = inputNodes;
    % compute x and y components of force
    Fx = Fmag * sin( thetaVIA( caseIDX ) ); % x component
    Fy = Fmag * cos( thetaVIA( caseIDX ) ); % y component
    % plug force values into array
    forces(:,2,caseIDX) = Fx; % x component
    forces(:,3,caseIDX) = Fy; % y component 
end

%% set target displacements

% displacement of each node
dTarget = zeros(nIO,2,nCases);
Target = zeros(nIO,2,nCases);

% create unscaled sinusoid with amplitude 1
unscaledSin = sin( linspace(0,2*pi,nIO*2+1) )



% populate array for each via point
for caseIDX = 1:nCases
    % set sinusoidal displacement for each output node
    for IOIDX = 1:nIO
        dTarget(IOIDX,1,caseIDX) = ampVIA(caseIDX) * unscaledSin(2*IOIDX);
    end
    Target(:,:,caseIDX) = coord_initial(outputNodes,[1 2]) + dTarget(:,:,caseIDX);
end











