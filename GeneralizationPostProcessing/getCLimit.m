function cMax = getCLimit(coordDeformed, conMatrix, lResting, elongationMax, forceIn)

% author: Reinier Kuppens and Ryan Lee

% this function evaluates the max scaling factor for force scaling

%% function inputs

% deformed coordinates
coordDeformed; % [Ncoord,3,Ncases] x,y,z coordinates of each node for each behavior
lResting; % [1] length of individual beam
elongationMax; % [1] max elongation in behavior
forceIn; % [1] force magnitude


%% function outputs

% aximum scale factor in case of force scaling
cMax = 0; % [1]

%% computation
coordLengths = coordDeformed(conMatrix(:,1),:,:)-coordDeformed(conMatrix(:,2),:,:);
linkLengths = sqrt(sum(coordLengths.^2,2));
elongation = (linkLengths - lResting);
maxE = max(abs(elongation), [], 'all');
cMax = forceIn*elongationMax/maxE;
end