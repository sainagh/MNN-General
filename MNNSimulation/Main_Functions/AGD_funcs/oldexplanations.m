% dLdu [NinputANDoutput*2,Ncases] stores derivatives of cost wrt output displacements (x and y) for each behavior
dLdu = zeros(NinputANDoutput*size(Target,2),Ncases);
% F [NDOFnodes*DOI,Ncases] stores force and moment inputs at each node, for each behavior
F = zeros(NDOFnodes*DOI,Ncases);
% u [NDOFnodes*DOI,Ncases] stores displacements (x and y) and rotations (in plane) for each node
u = zeros(NDOFnodes*DOI,Ncases);
% K [NDOFnodes*DOI,NDOFnodes*DOI] stiffness matrix
K = zeros(NDOFnodes*DOI,NDOFnodes*DOI);
% dFdu [NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases] derivative of force wrt displacements, defined as F and u above
dFdu = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases);
% dudF [NDOFnodes*DOI,NDOFnodes*DOI] derivative of displacement wrt force, dudF = -(dFdu)^-1, see derivation
dudF = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,Ncases);
% dFdK [NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI] derivative of force wrt stiffness, defined as F and K above; data structure is redundant, [NDOFnodes*DOI,NDOFnodes*DOI] is valid as there's zero terms for all times first and second indices are different, see derivation slides
dFdK = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
% dudK [NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI] derivative of displacement wrt stiffness, defined as u and K above
dudK = zeros(NDOFnodes*DOI,Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
% duoutdK [NinputANDoutput*2,Ncases,NDOFnodes*DOI,NDOFnodes*DOI] derivative of output displacements wrt stiffness, subset of dudK
duoutdK = zeros(NinputANDoutput*size(Target,2),Ncases,NDOFnodes*DOI,NDOFnodes*DOI);
% dLdK [NDOFnodes*DOI,NDOFnodes*DOI] derivative of cost function wrt stiffness, defined as L and K above
dLdK = zeros(NDOFnodes*DOI,NDOFnodes*DOI);
% dKdx [NDOFnodes*DOI,NDOFnodes*DOI,Nbeams] derivative of stiffness wrt control stiffness values
dKdx = zeros(NDOFnodes*DOI,NDOFnodes*DOI,Nbeams);
% dLdx [Nbeams] derivative of cost function wrt control stiffness values
dLdx = zeros(Nbeams,1);