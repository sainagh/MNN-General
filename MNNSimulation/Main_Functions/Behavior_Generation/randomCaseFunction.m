function [Target,forces, graphArrows, xyForces] = randomCaseFunction(LatticeGeometryStruct,BehaviorStruct)

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
elongationMax = BehaviorStruct.Elongation_max;
maxForce = BehaviorStruct.MaxForce;
threshold = BehaviorStruct.threshold;

%% computations


% obtain number of inputs and outputs
nIn = length(inputNodes);
nOut = length(outputNodes);

% 
Target = zeros(nOut, 2, nCases);
forces = zeros(nIn, 7, nCases);
for caseIndex = 1:nCases
  goodCase = 0;
  while goodCase == 0
    for node = 1:nOut
        nodeNumber = outputNodes(node);
        xDisp = ((rand*2)-1)*elongationMax; 
        yDisp = ((rand*2)-1)*elongationMax; 
        xPos = coord_initial(nodeNumber,1);
        yPos = coord_initial(nodeNumber,2);
        
        Target(node,:,caseIndex) = [xPos + xDisp, yPos + yDisp];
    end
    for node = 1:nIn
        nodeNumber = inputNodes(node);
        fA = (rand*2-1)*maxForce;
        fB = (rand*2-1)*maxForce;
        
        fx = cos(pi/4)*fA + sin(pi/4)*fB;
        fy = cos(pi/4)*fB - sin(pi/4)*fA;
        
        forces(node,[1,2,3],caseIndex) = [nodeNumber, fx, fy];
        xyForces(node, 2:4, caseIndex) = [nodeNumber, fx/maxForce, fy/maxForce];
    end
    if caseIndex ~= 1
        %subtract all the earlier force inputs from the new force input
        d = forces(:,[2,3],1:caseIndex-1)-forces(:,[2,3],caseIndex);
        mseF = sqrt( sum(d.^2, [1,2]) )/nIn; %calculate the mean squared error
        dx = Target(:,:,1:caseIndex-1)-Target(:,:,caseIndex);
        mseX = sqrt( sum(dx.^2, [1,2]))/nOut;
        if min(mseF) > threshold*maxForce && min(mseF) > threshold*elongationMax/sqrt(2)
            goodCase = 1;
        else
            caseIndex
%             plot(squeeze(mse));
        end
    else
        goodCase = 1;
    end
  end
    
end

for index = 1:size(forces,3)
    xforces1=find(forces(:,2,index));
    xforces1=[inputNodes' forces(xforces1,1,index) forces(xforces1,2,index)];
    yforces1=find(forces(:,3,index));
    yforces1=[inputNodes' forces(yforces1,1,index) forces(yforces1,3,index)];

    cellIndex = 2*index - 1;
    graphArrows{cellIndex} = xforces1;
    graphArrows{cellIndex + 1} = yforces1;
end
 
end
        
        
    