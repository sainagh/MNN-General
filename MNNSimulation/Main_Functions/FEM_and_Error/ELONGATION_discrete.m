function [c, ceq]=ELONGATION_discrete(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, discspaceidx,PossibleStiffnessArray)

% author: Pietro Sainaghi
% original program written by Reinier Kuppens and Ryan Lee

% this function evaluates the elongation of each link in a given deformed lattice, 
% and it computes the difference from the max and min elongations using a
% format that agrees with the matlab optimization toolbox formatting

%% function inputs

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases; % [1] number of force behaviors to optimize towards

% lattice structure - LatticeGeometryStruct
connectivity = LatticeGeometryStruct.connectivity; % [Nbeams,2] connectivity of nodes for each link

% design features structure - LinkPropertiesStruct
L1 = LinkPropertiesStruct.L1;  % [1] link length at rest
MaxLinkElongation = LinkPropertiesStruct.MaxLinkElongation; % [1] maximum deformation of link
MaxLinkCompression = LinkPropertiesStruct.MaxLinkCompression; % [1] maximum deformation of link
nonlinearStiffness = LinkPropertiesStruct.nonlinearStiffness; % [bool] wether to use nonlinearity in control stiffness

% link stiffness values
discspaceidx; % [Nbeams] stiffness value of each beam, selected from PossibleStiffnessArray
PossibleStiffnessArray;

%% function outputs

% excessElongation % [Nbeams, Ncases*2] difference (always negative) with respect to link length limits
% col 1: behavior 1 with respect to max elongation
% col 2: behavior 2 with respect to max elongation
% col 3: behavior 1 with respect to max compression
% col 4: behavior 2 with respect to max compression
          

%% Extract link stiffness values


for eachx = 1:length(discspaceidx)
    x(eachx) = PossibleStiffnessArray(discspaceidx(eachx));
end


%% compute deformed shape
if nonlinearStiffness
    coord_deformed=FEM_nonlinear(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
else
    coord_deformed=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
end
%% compute link lengths of deformed shapes

for behIDX = 1:Ncases
    for linkIDX = 1:size(connectivity,1)
        defomedLength(linkIDX,behIDX) = sqrt( (coord_deformed(connectivity(linkIDX,1),1)-coord_deformed(connectivity(linkIDX,2),1))^2 + (coord_deformed(connectivity(linkIDX,1),2)-coord_deformed(connectivity(linkIDX,2),2))^2 );
        excessElongation(linkIDX,behIDX) = defomedLength(linkIDX,behIDX) - (L1 + MaxLinkElongation);
        excessElongation(linkIDX,behIDX+Ncases) = -(defomedLength(linkIDX,behIDX) - (L1 - MaxLinkCompression));
    end
end

%% set outputs

ceq = [];
c = excessElongation;


end