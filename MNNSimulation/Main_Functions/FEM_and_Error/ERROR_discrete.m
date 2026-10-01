function [error]=ERROR_discrete(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,discspaceidx,PossibleStiffnessArray)

% author: Pietro Sainaghi
% original program written by Reinier Kuppens and Ryan Lee

% this function evaluates the mean squared error for a given lattice and a
% given set of behaviors

%% function inputs

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases; % [1] number of force behaviors to optimize towards
Target = BehaviorStruct.Target; % [NinputANDoutput, 2, Ncases] desired x and y coordinates for each output node
Elongation_max = BehaviorStruct.Elongation_max; % [1] maximum node displacement in output nodes

% lattice structure - LatticeGeometryStruct
connectivity = LatticeGeometryStruct.connectivity; % [Nbeams,2] connectivity of nodes for each link
Outputline = LatticeGeometryStruct.outputNodes; % [NinputANDoutput] indices of output nodes
NinputANDoutput = LatticeGeometryStruct.NinputANDoutput; % [1] number of input and output nodes in a given lattice
coord_initial = LatticeGeometryStruct.coord_initial; % [Ncoord,3] x,y,z coordinates of each node 

% design features structure - LinkPropertiesStruct
L1 = LinkPropertiesStruct.L1;  % [1] link length at rest
nonlinearStiffness = LinkPropertiesStruct.nonlinearStiffness; % [bool] wether to use nonlinearity in control stiffness

% finite truss structure - FEMStruct
F = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior

% optimization problem structure - OptimizerDataStruct
ForceScaling = OptimizerDataStruct.ForceScaling; % [boolean] yes or no to force scaling
doScaleMSE = OptimizerDataStruct.doScaleMSE; % [boolean] wether or not to do MSE scaling
MSEscale = OptimizerDataStruct.MSEscale;
MSEunits = OptimizerDataStruct.MSEunits;

    
% link stiffness values
discspaceidx; % [Nbeams] stiffness value of each beam, selected from PossibleStiffnessArray
PossibleStiffnessArray;
    
    
%% function output

error = 0; % MSE can be in length^2 or unitless depending on performance metrics



%% Assign x values based on indices
for eachx = 1:length(discspaceidx)
    x(eachx) = PossibleStiffnessArray(discspaceidx(eachx));
end

%% computations

if nonlinearStiffness
    coord_deformed=FEM_nonlinear(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
else
    coord_deformed=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
end

% %Peanalize links with deformatuons that are too large
% slope = 10000;
% dCoord = sum((coord_deformed - coord_initial).^3,2)-Elongation_max;
% dCoord(dCoord < 0) = 0;
% offset = slope*sum(dCoord);

%disp('ErrorRUN')
% boundsError = 0;
% boundsScaling = 100;
%disp('ErrorRUN')
%Calculate the scale factor of the force for optimal displacement
xa=(coord_deformed(Outputline,1:2,:) - coord_initial(Outputline,1:2));
xt = (Target - coord_initial(Outputline,1:2));
vectLength = numel(xt);
xar = zeros([vectLength,1]);
atr = zeros([vectLength,1]);
xar = reshape(xa,[vectLength,1]);
xtr = reshape(xt,[vectLength,1]);
c = sum(xar'*xtr)/sum(xar.^2);
cMax = getCLimit(coord_deformed,connectivity, L1, Elongation_max, 1);
if abs(c) > cMax
    c = sign(c)*cMax;
end
    
abserror=zeros(NinputANDoutput*Ncases,1);

% for i=1:Ncases  
%     abserror((i-1)*NinputANDoutput+(1:NinputANDoutput))=...
%         sqrt((sum((1000*coord_deformed(Outputline(1:NinputANDoutput),[1 2],i)-1000*Target(1:NinputANDoutput,[1 2],i)).^2,2)));
% end
for i=1:Ncases
    if ForceScaling == true
        abserror((i-1)*NinputANDoutput+(1:NinputANDoutput))=...
            ((sum((xa(:,:,i)*MSEunits*c-MSEunits*xt(:,:,i)).^2,2)));
    elseif ForceScaling == false
        abserror((i-1)*NinputANDoutput+(1:NinputANDoutput))=...
            ((sum((xa(:,:,i)*MSEunits-MSEunits*xt(:,:,i)).^2,2)));
    end
end
error=sum(abserror)/(NinputANDoutput*Ncases);
% error = sum((xa*1000000-xt*1000000).^2)/vectLength;
% error = sum((c*xa*1000-xt*1000).^2)/vectLength;

% convert to unitless if flag is active
if doScaleMSE
    error = error / MSEscale;
end



end