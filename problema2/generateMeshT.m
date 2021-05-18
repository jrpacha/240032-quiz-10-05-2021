function [nodes,elem] = generateMeshT(a,b,nDiv)
step=(b-a)/nDiv;
[X,Y] = meshgrid(a:step:b); % create a grid of points
    Z=zeros(size(X));
    [elemQ,nodes] = surf2patch(X,Y,Z);
    nodes=nodes(:,1:2); % to discard the third coordinate
% nodes=[X(:),Y(:)];
numElemQ=size(elemQ,1);
elem=zeros(2*numElemQ,3);
k=0;
for i=1:nDiv:numElemQ
    for j=i:1:i+nDiv-1
        k=k+1;
        elem(k,:)=[elemQ(j,4),elemQ(j,1),elemQ(j,2)];
        k=k+1;
        elem(k,:)=[elemQ(j,2),elemQ(j,3),elemQ(j,4)];
    end
end
