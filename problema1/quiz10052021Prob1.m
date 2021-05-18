%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quiz 10-05-2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 3;           %m
radius = 0.35;   %m   
Y = 5.0625e+06;  %N/m^2;
FA = 408.32;     %N;
FB = FA/2;
FC = FA/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%In this setting node A corresponds to node 8
%                node B corresponds to node 9
%                node C corresponds to node 7
%                node D corresponds to node 1
%                node E corresponds to node 5
numNodA = 8;
numNodB = 9;
numNodC = 7;
numNodD = 1;
numNodE = 5;

%In this setting elem AF corresponds to elem 14
numElemAF = 14;

%Geometry
nodes=[0,0;
       2,1;
       4,2;
       6,1;
       8,0;
       8,1;
       6,2;
       4,3;
       2,2;
       0,1];
   
nodes = nodes*h;   
       
elem=[1,2;
      2,3;
      3,4;
      4,5;
      5,6;
      6,7;
      7,8;
      8,9;
      9,10;
      10,1;
      10,2;
      2,9;
      9,3;
      3,8;
      3,7;
      7,4;
      4,6];
    
numNod=size(nodes,1);
numElem=size(elem,1);
dim=size(nodes,2);
  
numbering=1; %=0, neither labels nodes nor elements
plotElements(nodes, elem, numbering);

%Real constants: Materials and sections area
Area=pi*radius*radius;
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(dim*numNod,1);
Q=zeros(dim*numNod,1);
K=zeros(dim*numNod);

for e=1:numElem
    Ke=planeLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[2*elem(e,1)-1,2*elem(e,1),2*elem(e,2)-1,2*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads:
%node A (node 8):
nod = numNodA;
Q(dim*nod-1) = 0.0;
Q(dim*nod) = -FA; %N

%node B (node 9):
nod = numNodB;
Q(dim*nod-1)=0.0;
Q(dim*nod)= -FB;  %N

%node C (node 7):
nod = numNodC;
Q(dim*nod-1) = 0;
Q(dim*nod) = -FC; %N

%Boundary Conditions
fixedNods=[];
%node D (node 1):
nod = numNodD;
fixedNods = [fixedNods,dim*nod];     %(u1_y=0); 
u(dim*nod) = 0.0;

%node E (node 5):
nod = numNodE;
fixedNods= [fixedNods,dim*nod-1];   %(u5_x=0);
fixedNods = [fixedNods,dim*nod];     %(u5_y=0);
u(dim*nod-1) = 0.0;
u(dim*nod) = 0.0;

%Reduced system
freeNods=setdiff(1:dim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

ux = u(1:2:end);
uy = u(2:2:end);

%Post-process
%Show the original structure and the deformed one
%figure()
esc=5; %scale factor to magnify displacements
plotDeformedTruss(nodes,elem,u,esc);

%Reaction Forces:
RF = K*u-Q;

RFX = RF(1:2:end);
RFY = RF(2:2:end);

%Strain
L = nodes(elem(:,1),:)-nodes(elem(:,2),:);
LD = L+[ux(elem(:,1))-ux(elem(:,2)), uy(elem(:,1))-uy(elem(:,2))];
L = sqrt(L(:,1).^2+L(:,2).^2);  
LD = sqrt(LD(:,1).^2+LD(:,2).^2);
strain = abs(LD-L)./L;

%Print out displacements and reaction forces at nodes
%Fancy output: don't waste your time with this at the exams!
fprintf('\n%6s%10s%14s%14s%14s\n',...
    'NOD.','UX(m)','UY(m)','RX(N)','RY(N)')
fprintf('%4d%14.5e%14.5e%14.5e%14.5e\n',...
    [(1:numNod)', ux, uy, RFX, RFY]')

%Print out the the elements's strains
%Fancy output: don't waste your time with this at the exams!
fprintf('\n%6s%12s\n',...
    'ELEM.','STRAIN')
fprintf('%4d%16.5e\n',...
    [(1:numElem)', strain]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Print out the solution of part (a)
%Fancy output: don't waste your time with this at the exams!
fprintf('\n(a) Displacements of the point A (in m):\n')
fprintf('%6s%10s%14s\n', 'NOD.', 'UX(m)', 'UY(m)')
fprintf('%4d%14.5e%14.5e\n',numNodA, ux(numNodA), uy(numNodA))
fprintf('Hint. The x-displacement of  node D is: %.5e\n',ux(numNodD));

%Print out the solution of part (b)
%Fancy output: don't waste your time with this at the exams!
fprintf('\n(b) Reaction force at point D (in N):\n')
fprintf('%6s%10s%14s\n','NOD.','RX(N)','RY(N)')
fprintf('%4d%14.5e%14.5e\n',...
    [numNodD,RF(dim*numNodD-1),RF(dim*numNodD)]')

%Print out the solution of part (c)
%Fancy output: don't waste your time with this at the exams!
fprintf('\n(c) Maximal strain suffered by any element structure: %.5e\n',...
    max(strain))
fprintf('Hint. The strain for element AF is: %.5e\n',strain(numElemAF))