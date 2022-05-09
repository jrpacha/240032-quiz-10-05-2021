%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quiz 10-05-2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Question 2 (a, b)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 4;
nDiv = 20;
kc = 1.23;
TLR = 19;
q0 = -2;
p=[0.5,0.903];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nodes, elem] = generateMeshT(-L,L,nDiv);

numNodes= size(nodes,1);
numElem= size(elem,1);

numbering= 0; %= 1 shows nodes and element numbering
plotElements(nodes,elem,numbering)

%Select Boundary points
L1 = L-0.001;
L2 = -L+0.001;
indT= find(nodes(:,2) > L1); %indices of the nodes at the top boundary
indB= find(nodes(:,2) < L2); %indices of the nodes at the bottom boundary
indR= find(nodes(:,1) > L1); %indices of the nodes at the right edge
indL= find(nodes(:,1) < L2); %indices of the nodes at the left boundary

hold on %plot the points corresponding to the selected indexes
plot(nodes(indT,1),nodes(indT,2),'ok','lineWidth',1,'markerFaceColor',...
    'red','markerSize',5)
plot(nodes(indB,1),nodes(indB,2),'ok','lineWidth',1,'markerFaceColor',...
    'red','markerSize',5)
plot(nodes(indR,1),nodes(indR,2),'ok','lineWidth',1,'markerFaceColor',...
    'green','markerSize',5)
plot(nodes(indL,1),nodes(indL,2),'ok','lineWidth',1,'markerFaceColor',...
    'green','markerSize',5)
hold off

%Define the coefficients vector of the model equation
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=0;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
K=zeros(numNodes);    %global stiff matrix
F=zeros(numNodes,1);  %global internal forces vector
Q=zeros(numNodes,1);  %global secondary variables vector

for e = 1:numElem
    [Ke, Fe] = linearTriangElement(coeff,nodes,elem,e);
    rows= [elem(e,1); elem(e,2); elem(e,3)];
    cols= rows;
    K(rows,cols)= K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)= F(rows) + Fe;
    end
end

%Booundary Conditions
fixedNodes= [indL', indR'];         %fixed Nodes (global numbering)
freeNodes= setdiff(1:numNodes,fixedNodes); %free Nodes (global numbering)

%Natural B.C:
indBC = [indT', indB'];
Q=applyConstantNaturalBC(nodes,elem,indBC,q0,Q);

% Essential B.C.
u=zeros(numNodes,1);
u(indR)= TLR;
u(indL)= TLR;

%Reduced system
Fm = F(freeNodes) + Q(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);
Km = K(freeNodes,freeNodes);

%Compute the solution
um = Km\Fm;
u(freeNodes)= um;

%PostProcess: Compute secondary variables, table and plot results
Q = K*u - F;

% table = [(1:numNodes)',nodes(:,1),nodes(:,2),u,Q];
% fmt1 ='%4s%9s%14s%14s%14s\n';
% fmt2 ='%4d%14.5e%14.5e%14.5e%14.5e\n';
% fprintf(fmt1,'Node','X','Y','U','Q')
% fprintf(fmt2,table')

titol='Temperature Distribution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

for e=1:numElem
    vertexs= nodes(elem(e,:),:);
    [alphas,isInside] = baryCoord(vertexs,p);
    if (isInside >= 1)
        pElem = e;
        numNodElem= elem(e,:);
        tempP = alphas*u(numNodElem);
        break;
    end
end

%Print out solution of part (a)
%Fancy output: don't waste your time with this at the exams!
fprintf('\n')
fprintf('(a) Point P = (%f,%f) belongs to element #: %d\n',p,pElem)
fprintf(' Number of nodes of elem %d: %d, %d, %d\n',pElem,numNodElem)
fprintf(' Interpolated temperature at point P: %.5e%cC\n',tempP,char(176))
fprintf(' Hint. The minimum reached temperature is: %.5e%cC\n\n',...
    min(u),char(176));

%Find the nodes below the diagonal
indNodsBD = find(nodes(:,2) < nodes(:,1)+0.001); 
plotElements(nodes,elem,numbering)
hold on
plot(nodes(indNodsBD,1),nodes(indNodsBD,2),... % draw the nodes below diag.
    'ok','lineWidth',1,'markerFaceColor','red','markerSize',5)
hold off

%Compute the mean temperature for these nodes
avTempBD = sum(u(indNodsBD))/length(indNodsBD);

%Print out the solution of part (b)
%Fancy output: don't waste your time with this at the exams!
fprintf('(b) Mean of the temperatures of the nodes placed below \n')
fprintf('the diagonal, <T>_bd = %.5e%cC\n',avTempBD,char(176))