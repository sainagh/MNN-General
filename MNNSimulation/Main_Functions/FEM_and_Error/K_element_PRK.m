function [K_elem,axF]=K_element_PRK(LinkPropertiesStruct)
%% Construct stiffness matrix 2D
% axial direction = x
% perpendicular = y
% Moment = rotz

% Force acting axial
k44 = 1; % set to 1 because it's multiplied by the actual value in ERROR function

% Force acting perpendicular while rotation is fixed: utcross=0;
k55 = LinkPropertiesStruct.k55;

% Force acting perpendicular with rotation of tip allowed
k56 = LinkPropertiesStruct.k56;

% Moment acting on end while perpendicular translation of tip is allowed 
k66 = LinkPropertiesStruct.k66;

% Total stiffness matrix 
K_elem=[ k44    0     0 -k44    0     0;
           0  k55   k56    0 -k55  k56;
           0  k56   k66    0 -k56 k66/2;
        -k44    0     0  k44    0     0;
           0 -k55  -k56    0  k55  -k56;
           0  k56 k66/2    0 -k56   k66;];
       
axF=k44;
end