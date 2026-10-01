function BehaviorStruct = packageBehavior(LatticeGeometryStruct,BehaviorStruct)

% this function runs the function that generates a behavior input and
% target, and combines the outputs in a convenient structure


if BehaviorStruct.caseType == 1
    % generate sinusoid behavior
    [Target, ~, forces, ~, graphArrows, xyForces] = pushShearSinWaveFunction...
        (LatticeGeometryStruct,BehaviorStruct);
    
elseif BehaviorStruct.caseType == 2
    % generate randomized behavior
    
    % type of rng
    switch BehaviorStruct.RNGtype
        case 'dateandtime'
            t= datetime;
            tt = datenum(t);
            rng(tt);
        case 'deterministic'
            rng(BehaviorStruct.RNG_seed_det);
        case 'nocontrol'
            % do nothing
    end

    [Target, forces, graphArrows, xyForces] = randomCaseFunction...
        (LatticeGeometryStruct,BehaviorStruct);
    
    elseif BehaviorStruct.caseType == 4
        % generate generalizable behavior
        [Target, forces] = generalizedCaseFunction(LatticeGeometryStruct,BehaviorStruct);
        
end

% store behavior data into structure BehaviorStruct
% target node coordinates
BehaviorStruct.Target = Target;
% forces locations in DOF vector format
BehaviorStruct.forces = forces;