%................................................................

% MATLAB codes for Finite Element Analysis
% problem17_planestrain.m
% 2D problem: thin plate in tension
% original create by antonio ferreira 2008
% modified by J.Tu   2020

% clear memory
clear all;colordef white;clf

% materials
E  = 10e7;     poisson = 0.30;  

% matriz C
C=E/((1+poisson)*(1-2*poisson))*[(1-poisson) poisson 0;poisson (1-poisson) 0;0 0 ((1-2*poisson)/2)];

% load
P = 1e6;

%Mesh generation
Lx=5;
Ly=1;
numberElementsX=10;
numberElementsY=5;
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(Lx,Ly,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% calculation of the system stiffness matrix
stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C,1,1);

% boundary conditions 
fixedNodeX=find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY=find(nodeCoordinates(:,2)==0);  % fixed in YY
prescribedDof=[fixedNodeX; fixedNodeY+numberNodes];

% force vector (distributed load applied at xx=Lx)
force=zeros(GDof,1);
rightBord=find(nodeCoordinates(:,1)==Lx);
force(rightBord)=P*Ly/numberElementsY;
force(rightBord(1))=P*Ly/numberElementsY/2;
force(rightBord(end))=P*Ly/numberElementsY/2;

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
scaleFactor=10;

% deformed shape
figure
drawingField(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,'Q4',UY);%U XX
hold on
drawingMesh(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,'Q4','k-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
colorbar
title('U XX  (on deformed shape)')
axis off

% stresses at nodes
stresses2D(GDof,numberElements,elementNodes,numberNodes,...
    nodeCoordinates,displacements,UX,UY,C,scaleFactor)