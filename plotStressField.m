function plotStressField(numberElements,numberNodes,...
  elementNodes,stress_NP,nodeCoordinates,UX,UY,...
  scaleFactor)

% stress_node_global = zeros(4,numberElements)
global stress_node_average
global stress_Node_Total
for h = 1:numberNodes

% find this Node in the structure.
    [i,j]=find(elementNodes==h)

% this node is meshed by how many elements.	
	node_shared = size(i,1)
	
% node infortaion: which include this node in which element,
% and node NO. in the element.
	node_Infor = [i,j]
	
	
	stress_Node_Total = zeros(4,1)
	for e = 1:node_shared
	
    %determine the information of Node (in which element, and the NO.)	
	node_Localation = node_Infor(e,:)

    %export the node stress of each element
	stress_Node_Element(:,1) = stress_NP(:,node_Localation(2),node_Localation(1))
	
	%add all node stress from related element.
    stress_Node_Total(:,1)   = stress_Node_Total(:,1) + stress_Node_Element(:,1)
	end 

	% average the node stress from related element
	stress_node_average(1,h) = (sum(stress_Node_Total(1,1)))/node_shared
	stress_node_average(2,h) = (sum(stress_Node_Total(2,1)))/node_shared
	stress_node_average(3,h) = (sum(stress_Node_Total(3,1)))/node_shared
	stress_node_average(4,h) = (sum(stress_Node_Total(4,1)))/node_shared	
end

% To imput to function 'drawingField', transpose the matrix
stress_node_average_T_X  = (stress_node_average(1,:))'
stress_node_average_T_Y  = (stress_node_average(2,:))'
stress_node_average_T_Z  = (stress_node_average(4,:))'
stress_node_average_T_XY = (stress_node_average(3,:))'

for h = 1:numberNodes

	stress_node_average_T_VMS(h) = ...
	((((stress_node_average(1,h)-stress_node_average(2,h))^2)+ ...
	((stress_node_average(2,h)-stress_node_average(4,h))^2)+ ...
	((stress_node_average(4,h)-stress_node_average(1,h))^2)+ ...
	6*(stress_node_average(3,h)^2))/2)^(0.5)

end

% deformed shape
figure
plotField2D(nodeCoordinates+scaleFactor*[UX UY], ...
elementNodes,stress_node_average_T_X,numberElements)%Stress XX
hold on
%plotMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,1)
drawingMesh(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,'Q4','k-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
colorbar
title('Sigma radial direction stress (on deformed shape)')
axis off

% deformed shape
figure
plotField2D(nodeCoordinates+scaleFactor*[UX UY], ...
elementNodes,stress_node_average_T_Y,numberElements)%Stress YY
hold on
%plotMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,1)
drawingMesh(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,'Q4','k-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
colorbar
title('Sigma axial direction stress (on deformed shape)')
axis off

% deformed shape
figure
plotField2D(nodeCoordinates+scaleFactor*[UX UY], ...
elementNodes,stress_node_average_T_XY,numberElements)%Stress XY
hold on
%plotMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,1)
drawingMesh(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,'Q4','k-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
colorbar
title('Sigma XY stress (on deformed shape)')
axis off

% deformed shape
figure
plotField2D(nodeCoordinates+scaleFactor*[UX UY], ...
elementNodes,stress_node_average_T_VMS,numberElements)%Stress Von Mises Stress
hold on
%plotMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,1)
drawingMesh(nodeCoordinates+scaleFactor*[UX UY],...
    elementNodes,'Q4','k-');
drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
colorbar
title('Von Mises stress (on deformed shape)')
axis off