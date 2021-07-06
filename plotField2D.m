% plot field contour

function plotField2D(nodeCoordinates, ...
elementNodes,value,numberElements)

%---------------------------------------------
hold on;  
for e=1:numberElements  
    i=elementNodes(e,:);%Connectivity Matrix  
    x=nodeCoordinates(i,1);  
    y=nodeCoordinates(i,2);  
    s=(value(i))';  
    fill(x,y,s,'FaceColor','interp');  
end  

shading interp;
colorbar 
axis equal