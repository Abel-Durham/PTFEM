function [vonMisesStress]=vonMisesStress(StressMatrix)

% this function is used to calculate the von Mises Stress.


% von Misesstresses 
vonMisesStress = zeros(size(StressMatrix,2),size(StressMatrix,3));
numberPoints   = size(StressMatrix,2)
numberElements = size(StressMatrix,3)
  
  
  
for e=1:numberElements                           

	for q=1:numberPoints
		
		stress = StressMatrix(:,q,e)

		equivalentStress = ((((stress(1)-stress(2))^2)+ ...
							 ((stress(2)-stress(4))^2)+ ...
							 ((stress(4)-stress(1))^2)+ ...
							 6*(stress(3)^2))/2)^(0.5)
			
		vonMisesStress(q,e) = equivalentStress
	end
end	