%................................................................

 function [stress_GP,strain_GP]=stresses2D_GP(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,...
    displacements,UX,UY,C,scaleFactor)

% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');

% stresses at Gaussian Point
stress_GP=zeros(4,size(elementNodes,2),numberElements);
  
for e=1:numberElements                           
  indice=elementNodes(e,:); 
  elementDof=[ indice indice+numberNodes ];   
  nn=length(indice);
  for q=1:size(gaussWeights,1)                        
    pt=gaussLocations(q,:);                            
    wt=gaussWeights(q);                              
    xi=pt(1);
    eta=pt(2);
% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

% the global coordinates of gaussion point
	gc=shapeFunction'*(nodeCoordinates(indice,:))

%  B matrix
    B=zeros(4,2*nn);
    B(1,1:nn)       = XYderivatives(:,1)';
    B(2,nn+1:2*nn)  = XYderivatives(:,2)';
    B(3,1:nn)       = XYderivatives(:,2)';
    B(3,nn+1:2*nn)  = XYderivatives(:,1)';
	B(4,1:nn) = (shapeFunction/gc(1))'
    
% element deformation 
    strain=B*displacements(elementDof);
	strain_GP(:,q,e)=strain;
    stress_GP(:,q,e)=C*strain;
	end
end	