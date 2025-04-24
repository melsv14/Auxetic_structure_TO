%%%% 3D Iso-XFEM topology optimization code by Meisam Abdi,University of Nottingham %%%%
% This code generates topology optimization solutions for 3D structures
% using Iso-XFEM method.
function BESO2D(L,h,nelx,nely,volfrac,E,nu,er,rmin)
% L: length, W: width, h: height
% nelx: number of elements in x direction
% nely: number of elements in y direction
% nelz: number of elements in z direction
% er: volume evolution rate
% volfrac: volume constraint
% E: Young's modulus, nu: Poisson’s ratio
% f: magnitude of load
% BESO3D(40,2,20,20,2,10,0.5,1,0.3,0.01,1.5) % example
[con,coord]=concoord(nelx,nely,L,h); % Element connectivity and Nodes' coordinate matrix
nel=size(con,1);nnd=size(coord,1); %number of nodes and elements
% Element elasticity matrix
D=E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; % Element elasticity matrix
% initial design
x=ones(nel,1); % set all design variables as 1 (full solid design domain)
% Display and save the start mesh (full solid design domain)
figure (1)
for i=1:nel
for j=1:4
xcoord(j,i)=coord(con(i,j),1);
ycoord(j,i)=coord(con(i,j),2);
end
end
axis equal; 
axis([min(coord(:,1)) max(coord(:,1)) min(coord(:,2)) max(coord(:,2)) ]);
axis off;
h=patch(xcoord,ycoord,'red');
set(h,'LineWidth',1.5)
set(h,'EdgeColor','red');
% set(h,'FceColor',[.85 .85 .85]);
hgsave('0')
% Force and Boundary Conditions
U=zeros(2*nnd,1);
F=zeros(2*nnd,1);
bottomend=2*(nelx+1)*(nely+1)-2*nely;
% fixeddofs=[1 2 bottomend];
fixeddofs=[1:2*(nely+1)];
alldofs = [1:2*nnd];
freedofs = setdiff(alldofs,fixeddofs);
%bottomend should be 2*(nelx+1)*(nely+1)-2*nely... without -2*nely is for
%upperend....!!!!!
% bottomend=2*(nelx+1)*(nely+1)-2*nely; %DoF number of the bottom of the free end
middleend=2*(nelx+1)*(nely+1) -nely; %DoF number of the middle of the free end
fdof= middleend;
% fdof=(nelx+1)*(nely+1)+(nely+1);
% fdof=(nelx+1)*(nely+1)+1;
F(fdof)=-1;
% elements connected to each node
ncon=zeros(nnd,4);
for i=1:nnd
row=[];col=[];
[row,col] = find(con==i);
d=size(row,1);
ncon(i,1:d)=row;
end
%volume and stiffness of the elements
for el=1:nel
elcoord=coord(con(el,:)',:);
[~,Ae1(el,1)] = convhulln([elcoord(1,:);elcoord(2,:);elcoord(3,:)]);
Ke1(:,:,el)=stiffnessmat(D,el,coord,con); % element stiffness at start mesh
end
A=sum(sum(Ae1)); % area of the whole structure
Ke=Ke1;Ae=Ae1; % set the current stifnees and volume as those of the initial design
if rmin>=0
[nbnodes,nbdis]=nodedis(coord,rmin);
end
volfracit=1; % defining initial volume fraction of the design
it=0;
results=[];
resultsIt=[];
resultsSE=[];
resultsV=[];
while (((volfracit-volfrac)/volfrac)>=0.001)
    it=it+1;
if it>1; olddc=dc; oldndc=ndc; end
volfracit=max(volfracit*(1-er),volfrac); % calculate volume of the next design
[K]=stiffness(nel,nnd,con,Ke,x); % Global stiffness matrix
U(freedofs,1)=K(freedofs,freedofs)\F(freedofs,1); % FEA for initial structure
for el = 1:nel
edof = [2*con(el,1)-1; 2*con(el,1); 2*con(el,2)-1;
2*con(el,2); 2*con(el,3)-1; 2*con(el,3); 2*con(el,4)-1;
2*con(el,4)];
Ue=U(edof,1);
ese(el,1)=0.5*Ue'*Ke(:,:,el)*x(el)*Ue; % element strain energy
dc(el,1) = ese(el,1)/Ae1(el,1); % element strain energy density (elemental sensiytivities)
end
SE=sum(sum(ese));
for i=1:nnd
a=find(ncon(i,:)>0);
b=dc(ncon(i,a),1);
% ndc(i,1)=mean(b);
ndc(i,1)=sum(b)/numel(a); % nodal sensitivities
end
% Filtering The Sensitivities
if rmin>0
[ndc,dc]=filterBESO(nbnodes,nbdis,ndc,dc,rmin,con);
end
if it>1; dc=(dc+olddc)/2; ndc=(ndc+oldndc)/2; end % stabilization of evolutionary process
vf(it)=sum(sum(Ae))/A; % volume fraction of current design
% BESO Algorithm
[x]=ADDDEL(volfracit,A,dc,x,Ae);
sol=x==1;
vol(it)=sum(Ae1(sol,1));
%%%%%%%%%%%%%%%%%%%%%%%%% visualization of the topology
surfcon=[];xcoord=[];ycoord=[];
x2coord=[];
y2coord=[];
% for i=1:nel
% if x(i)==1
% surfcon(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5])
% con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
% end
% end
figure (it)
m=0;
z=0;
for i=1:nel
    if x(i)==1
        m=m+1;
        for j=1:4
            xcoord(j,m)=coord(con(i,j),1);
            ycoord(j,m)=coord(con(i,j),2);
        end
    else
        z=z+1;
        for j=1:4
            x2coord(j,z)=coord(con(i,j),1);
            y2coord(j,z)=coord(con(i,j),2);
        end
    end
end
% figure (it)
h2=patch(xcoord,ycoord,'red');
h3=patch(x2coord,y2coord,'blue');
axis equal; axis([min(coord(:,1)) max(coord(:,1)) min(coord(:,2)) max(coord(:,2))]);axis off;
set(h2,'LineWidth',1.5)
set(h2,'EdgeColor','red');
set(h3,'LineWidth',1.5)
set(h3,'EdgeColor','blue');
% camlight
% % colormap(gray(100))
% lighting gouraud
% hgsave(num2str(it))
if it<200
close(it);
end
disp([' It.: ' sprintf('%4i',it) ' Obj.: ' sprintf('%10.4f',SE) 'Volfrac.: ' sprintf('%6.3f',vol(it)/A )]);
resultsIt(it)=[it];
resultsSE(it)=[SE];
resultsV(it)=[vol(it)/A ];
end
for w=1:10
 it=it+1;
if it>1; olddc=dc; oldndc=ndc; end
volfracit=max(volfracit*(1-er),volfrac); % calculate volume of the next design
[K]=stiffness(nel,nnd,con,Ke,x); % Global stiffness matrix
U(freedofs,1)=K(freedofs,freedofs)\F(freedofs,1); % FEA for initial structure
for el = 1:nel
edof = [2*con(el,1)-1; 2*con(el,1); 2*con(el,2)-1;
2*con(el,2); 2*con(el,3)-1; 2*con(el,3); 2*con(el,4)-1;
2*con(el,4)];
Ue=U(edof,1);
ese(el,1)=0.5*Ue'*Ke(:,:,el)*x(el)*Ue; % element strain energy
dc(el,1) = ese(el,1)/Ae1(el,1); % element strain energy density (elemental sensiytivities)
end
SE=sum(sum(ese));
for i=1:nnd
a=find(ncon(i,:)>0);
b=dc(ncon(i,a),1);
% ndc(i,1)=mean(b);
ndc(i,1)=sum(b)/numel(a); % nodal sensitivities
end
% Filtering The Sensitivities
if rmin>0
[ndc,dc]=filterBESO(nbnodes,nbdis,ndc,dc,rmin,con);
end
if it>1; dc=(dc+olddc)/2; ndc=(ndc+oldndc)/2; end % stabilization of evolutionary process
vf(it)=sum(sum(Ae))/A; % volume fraction of current design
% BESO Algorithm
[x]=ADDDEL(volfracit,A,dc,x,Ae);
sol=x==1;
vol(it)=sum(Ae1(sol,1));
%%%%%%%%%%%%%%%%%%%%%%%%% visualization of the topology
surfcon=[];xcoord=[];ycoord=[];
x2coord=[];
y2coord=[];
% for i=1:nel
% if x(i)==1
% surfcon(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5])
% con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
% end
% end
figure (it)
m=0;
z=0;
for i=1:nel
    if x(i)==1
        m=m+1;
        for j=1:4
            xcoord(j,m)=coord(con(i,j),1);
            ycoord(j,m)=coord(con(i,j),2);
        end
    else
        z=z+1;
        for j=1:4
            x2coord(j,z)=coord(con(i,j),1);
            y2coord(j,z)=coord(con(i,j),2);
        end
    end
end
h2=patch(xcoord,ycoord,[0.6350 0.0780 0.1840]);
h3=patch(x2coord,y2coord,[0 0.4470 0.7410]);
axis equal; axis([min(coord(:,1)) max(coord(:,1)) min(coord(:,2)) max(coord(:,2))]);axis off;
set(h2,'LineWidth',1.5)
% set(h2,'FaceColor',[.85 .85 .85])
set(h2,'EdgeColor','red');
h2.EdgeColor = [0.6350 0.0780 0.1840];
set(h3,'LineWidth',1.5)
% set(h2,'FaceColor',[.85 .85 .85])
% set(h3,'EdgeColor','blue');
h3.EdgeColor=[0 0.4470 0.7410];
% camlight
% colormap(gray(100))
% lighting gouraud
hgsave(num2str(it))
disp([' It.: ' sprintf('%4i',it) ' Obj.: ' sprintf('%10.4f',SE) 'Volfrac.: ' sprintf('%6.3f',vol(it)/A )]);   
resultsIt(it)=[it];
resultsSE(it)=[SE];
resultsV(it)=[vol(it)/A ];
end
results(:,1)=resultsIt;
results(end+1:end+it,1)=resultsSE;
results(end+1:end+it,1)=resultsV;
% dlmwrite('resultsBridge220110Ratio10WithFilterEd2.txt', results, 'precision', 15);
end