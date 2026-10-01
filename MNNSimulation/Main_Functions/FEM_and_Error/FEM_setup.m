function [DOFFinal,Final, F, DOFnodes,k_base,R6,RT6,DOF,Degrees_per_element] = FEM_setup(LinkPropertiesStruct,LatticeGeometryStruct,BehaviorStruct)


% author: Pietro Sainaghi
% original program written by Reinier Kuppens and Ryan Lee


%% function inputs

% design features structure - LinkPropertiesStruct
L1 = LinkPropertiesStruct.L1;  % [1] link length at rest

% lattice structure - LatticeGeometryStruct
Nnodes = LatticeGeometryStruct.Ncoord;
DOI = LatticeGeometryStruct.DOI;
NinputANDoutput = LatticeGeometryStruct.NinputANDoutput;
Outputline = LatticeGeometryStruct.outputNodes;
coord_initial = LatticeGeometryStruct.coord_initial;
bound = LatticeGeometryStruct.bound;
connectivity = LatticeGeometryStruct.connectivity;
Nbeams = LatticeGeometryStruct.Nbeams;

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases;
Target = BehaviorStruct.Target;
forces = BehaviorStruct.forces;

%% finite element setup

% degrees of freedom
DOF=DOI*Nnodes; %Total number of system DOF

% force vector
F=zeros(DOF,Ncases);

% store error
Error=zeros(3,1); %Error xyz
abserror=zeros(NinputANDoutput*Ncases,1); %normalized error for each output and case


%This sets the initial eror for each case
if Ncases == 1
    abserror(1:NinputANDoutput)=sqrt(sum((coord_initial(Outputline(1:NinputANDoutput),[1 2],1)-Target(1:NinputANDoutput,[1 2],1)).^2,2));
else
    for i=1:Ncases
        abserror((i-1)*NinputANDoutput+(1:NinputANDoutput))=sqrt(sum((coord_initial(Outputline(1:NinputANDoutput),[1 2],1)-Target(1:NinputANDoutput,[1 2],i)).^2,2));
    end
end

%Error(1) holds the average of the mean square error before
%optimization
Error(1)=sum(abserror.^2)/(NinputANDoutput*2);
% Make Force Vector
for i=1:size(forces,1)
    for j=1:Ncases
        %Holds the
        startIndex = (forces(i,1)-1)*DOI +1; %index from 1
        endIndex   =  startIndex + (DOI - 1); %To correct for index from 1
        F( [startIndex:endIndex], j)=forces(i,2:DOI+1,j);
    end
end

% fixed elements
Nfixed=length(bound);
prescribedDof=zeros(Nfixed,DOI);
for i=1:Nfixed
    prescribedDof(i,:)=(bound(i)-1)*DOI+1:(bound(i)-1)*DOI+DOI;
end
%X=coord_initial(:,1);Y=coord_initial(:,2);Z=coord_initial(:,3);
Final=setdiff(transpose(1:DOF),prescribedDof);
DOFFinal=length(Final);

%% fixed values

% undisturbed positions
pos1=coord_initial(connectivity(:,1),:,1); %xyz pos of link startPt
pos2=coord_initial(connectivity(:,2),:,1); %xyz pos of link endPt

% degrees of freedom in each element
Degrees_per_element=[DOI*connectivity(:,1)-2 DOI*connectivity(:,1)-1 DOI*connectivity(:,1) DOI*connectivity(:,2)-2 DOI*connectivity(:,2)-1 DOI*connectivity(:,2)] ;

% geometric parameters
RX = (pos2(:,1)-pos1(:,1))/L1;
RYX = (pos2(:,2)-pos1(:,2))/L1;
RZX = (pos2(:,3)-pos1(:,3))/L1;
D = sqrt(RX.*RX + RYX.*RYX);
RXy = -RYX./D;
RY = RX./D;
RZY = zeros(Nbeams,1);
RXz = -RX.*RZX./D;
RYz = -RYX.*RZX./D;
RZ = D;
mat2=zeros(3,3,Nbeams);
R12=zeros(12,12,Nbeams);

% material matrices
for i=1:Nbeams
    mat2(:,:,i) = [RX(i) RYX(i) RZX(i) ;RXy(i) RY(i) RZY(i) ;RXz(i) RYz(i) RZ(i)];
    R12(:,:,i) = [mat2(:,:,i) zeros(3,9); zeros(3) mat2(:,:,i) zeros(3,6);
        zeros(3,6) mat2(:,:,i) zeros(3);zeros(3,9) mat2(:,:,i)];
end
Relevant_R=[1 2 6 7 8 12];%Pulls off values that participate
R6=R12(Relevant_R, Relevant_R,:);
RT6=permute(R6,[2 1 3]);         %Rotation Matrix 3D
DOFnodes=setdiff(1:Nnodes,bound);
k_base=zeros(6,6,Nbeams);
for i=1:Nbeams
    [k_base(:,:,i),~]=K_element_PRK(LinkPropertiesStruct);    %Stiffness matrix?
end
[~,Kaxial]=K_element(LinkPropertiesStruct);

