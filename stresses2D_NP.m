 function [stress_NP]=stresses2D_NP(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,...
    displacements,UX,UY,C,scaleFactor)

stress = stresses2D_GP(GDof,numberElements,elementNodes,numberNodes,...
    nodeCoordinates,displacements,UX,UY,C,scaleFactor)


% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');

% local coordinates of each Gaussian Point
xg1=gaussLocations(1,1)
yg1=gaussLocations(1,2)
xg2=gaussLocations(2,1)
yg2=gaussLocations(2,2)
xg3=gaussLocations(3,1)
yg3=gaussLocations(3,2)
xg4=gaussLocations(4,1)
yg4=gaussLocations(4,2)

% matrxi for calculate the coefficient matrix 
% at Gaussian Point
mg = [xg1*yg1,xg1,yg1,1;xg2*yg2,xg2,yg2,1;xg3*yg3,xg3,yg3,1;xg4*yg4,xg4,yg4,1;]

% local coordinates of each Node Point
stressPoints=[-1 -1;1 -1;1 1;-1 1];
xn1=stressPoints(1,1)
yn1=stressPoints(1,2)
xn2=stressPoints(2,1)
yn2=stressPoints(2,2)
xn3=stressPoints(3,1)
yn3=stressPoints(3,2)
xn4=stressPoints(4,1)
yn4=stressPoints(4,2)

% matrxi for calculate the coefficient matrix 
% at Gaussian Point
mn = [xn1*yn1,xn1,yn1,1;xn2*yn2,xn2,yn2,1;xn3*yn3,xn3,yn3,1;xn4*yn4,xn4,yn4,1;]


stress_NP = zeros(4,size(elementNodes,2),numberElements);
  
for e=1:numberElements                           

  for q=1:size(gaussWeights,1)
%  Extract the stress component of each element 
%  at the Gauss integration point      
      stressXX(q) = stress(1,q,e)
	  stressYY(q) = stress(2,q,e)
	  stressXY(q) = stress(3,q,e)
	  stressZZ(q) = stress(4,q,e)	  
  end
  
% Calculate the matrix of interpolation coefficients
% for each stress component

	  cXX = inv(mg)*stressXX(:)
	  cYY = inv(mg)*stressYY(:)
	  cXY = inv(mg)*stressXY(:)
	  cZZ = inv(mg)*stressZZ(:)
	  
	  
%  calculate the stress component of each element 
%  at the node point
	  stressXX_N = mn*cXX
	  stressYY_N = mn*cYY	  
	  stressXY_N = mn*cXY
	  stressZZ_N = mn*cZZ
  
  for q=1:size(elementNodes,2) 	 
      stress_NP(1,q,e)= stressXX_N(q)   
      stress_NP(2,q,e)= stressYY_N(q) 	 
      stress_NP(3,q,e)= stressXY_N(q) 	 
      stress_NP(4,q,e)= stressZZ_N(q) 
	 	  
  end	 	  
end