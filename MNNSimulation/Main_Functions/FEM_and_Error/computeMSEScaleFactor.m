function [MSEscale] = computeMSEScaleFactor(LatticeGeometryStruct,BehaviorStruct)


%% function inputs

% lattice structure - LatticeGeometryStruct
coord_initial = LatticeGeometryStruct.coord_initial;
Outputline = LatticeGeometryStruct.outputNodes;
NinputANDoutput = LatticeGeometryStruct.NinputANDoutput;

% behavior structure - BehaviorStruct
Target = BehaviorStruct.Target;
Ncases = BehaviorStruct.Ncases;

%% computations
% compute undeformed error in each coordinate
undeformedError = (Target - coord_initial(Outputline,1:2));

% compute square of coordinates for each output point
for casenum = 1:Ncases
    absUndeformedError( (casenum-1) * NinputANDoutput + (1:NinputANDoutput) ) = sum( (1000000*undeformedError(:,:,casenum)).^2,2 );
end

% compute MSE scale factor
MSEscale = sum(absUndeformedError)/(NinputANDoutput*Ncases);