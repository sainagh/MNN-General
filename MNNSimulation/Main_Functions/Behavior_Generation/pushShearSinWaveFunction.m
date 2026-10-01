function [Target,Outputline,forces, dx, graphArrows, xyForces] = pushShearSinWaveFunction(LatticeGeometryStruct,BehaviorStruct)

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
maxDLength = BehaviorStruct.dx;
maxForce = BehaviorStruct.MaxForce;

%% initial formatting
    dx=maxDLength;%*Nlayers;
    nIO = length(inputNodes);
    forces=zeros(nIO,7,nCases);
    %% Input
        
        forces(:,1,1)= inputNodes;
        forces(:,1,2)= inputNodes;
        %Case 1: Compression in positive x direction
            Fx=maxForce/nIO; %Newton
            forces(:,2,1)=Fx*ones(nIO,1);
        %Case 2: Shear in positive y direction
            Fy=maxForce/nIO; %Newton
            forces(:,3,2)=Fy*ones(nIO,1);
    
       for index = 1:size(forces,3)
            xforces1=find(forces(:,2,index));
            xforces1=[xforces1 forces(xforces1,1,index) forces(xforces1,2,index)];
            yforces1=find(forces(:,3,index));
            yforces1=[yforces1 forces(yforces1,1,index) forces(yforces1,3,index)];
            
            cellIndex = 2*index - 1;
            graphArrows{cellIndex} = xforces1;
            graphArrows{cellIndex + 1} = yforces1;
       end
       xyForces = zeros(nIO,4,nCases);
       xyForces(:,2) = inputNodes;
       for index = 1:nCases
           xyForces(:,2,index) = inputNodes;
           xyForces(:, 3:4, index) = forces(:,2:3,index);
       end          
           
    %% Output
        dx=maxDLength/2;%*Nlayers;
        dTarget=zeros(nIO,DOI-1,nCases);
        for i=1:nIO
            rad=(i*2-1)/(2*nIO)*(2*pi);
            dTarget(i,1,1)=sin(rad)*dx;
            dTarget(i,1,2)=-sin(rad)*dx;
        end

        Outputline= outputNodes; %Which row(s) of the coord matrax are the outputs

        Target=zeros(nIO,DOI-1,nCases);
        Target(:,:,1)=coord_initial(Outputline,[1 2])+dTarget(:,:,1); %Target is position of nodes + displacement in target Pos 2
        Target(:,:,2)=coord_initial(Outputline,[1 2])+dTarget(:,:,2); %Target is position of nodes + displacement in target pos 2

        
%% reformat bumps

%Swap the targets so that push in X leads to Top bump OUT
Target = Target(:,:,nCases:-1:1);
% Target(:,1,:) = Target(:,1,:)+1e-6;
forces = forces(:,:,1:nCases);
    %% returnSection
%     Target
%     Outputline
%     forces
%     dx
    %These are used inplot function
end