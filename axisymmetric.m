%................................................................

% MATLAB codes for Finite Element Analysis
% problem17_axisymmetric.m
% 2D aixsymmetric problem
% original create by antonio ferreira 2008
% clear memory
clc
clear all;colordef white;clf

% materials
E  = 200e6;     poisson = 0.29;  

% matriz C
C=(E*(1-poisson))/((1+poisson)*(1-2*poisson))*[1 poisson/(1-poisson) 0 poisson/(1-poisson); ...
poisson/(1-poisson) 1 0 poisson/(1-poisson); 0 0 ((1-2*poisson)/(2*(1-poisson))) 0; ... 
poisson/(1-poisson) poisson/(1-poisson) 0 1]

% load
P = 5e3;

%FE model information
nodeCoordinates = nCoordinates;
elementNodes = nodeTopology;
nodeCoordinates = nodeCoordinates(1:end,2:3);
elementNodes = elementNodes(1:end,2:5);
numberElements = size(elementNodes(:,1));
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% calculation of the system stiffness matrix
stiffness=formStiffness2D_axisymmetric(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C);
	
% boundary conditions 
fixedNodeX=find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY=find(nodeCoordinates(:,2)==0);  % fixed in YY
prescribedDof=[fixedNodeX; fixedNodeY+numberNodes];

% force vector 
force=zeros(GDof,1);
loadingIndex = index;
force(int32(loadingIndex(1:end,2))) = loadingIndex(1:end,1);
force = force*P;



% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% displacements
disp('Displacements')
jj=1:GDof; format
f=[jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)
UX=displacements(1:numberNodes);
UY=displacements(numberNodes+1:GDof);
scaleFactor=1;

% deformed shape
figure
plotField2D(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,UX,numberElements);%U XX
hold on
drawingMesh(nodeCoordinates+scaleFactor*[UX UY],...
   elementNodes,'Q4','k-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
colorbar
title('U XX  (on deformed shape)')
axis off

% stress & strain at Gaussian Point
[stress_GP,strain_GP] = stresses2D_GP(GDof,numberElements,elementNodes,numberNodes,...
    nodeCoordinates,displacements,UX,UY,C,scaleFactor)
	
% stress at Node Point
stress_NP = stresses2D_NP(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,...
    displacements,UX,UY,C,scaleFactor)
	
% plot the stress field
plotStressField(numberElements,numberNodes,...
  elementNodes,stress_NP,nodeCoordinates,UX,UY,...
  scaleFactor)
 
% calculation of EQUIVALENT STRESS
stress_GP_vonMises=vonMisesStress(stress_GP)
stress_NP_vonMises=vonMisesStress(stress_NP)

%max von Mises Stress
maxVonMisesStress = max(max(stress_GP_vonMises))
