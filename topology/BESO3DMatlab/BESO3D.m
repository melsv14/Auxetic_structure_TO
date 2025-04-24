%%%% 3D Iso-XFEM topology optimization code by Meisam Abdi,University of Nottingham %%%%
% This code generates topology optimization solutions for 3D structures
% using Iso-XFEM method.
function BESO3D(L,W,h,nelx,nely,nelz,volfrac,E,NU,er,rmin)
% L: length, W: width, h: height
% nelx: number of elements in x direction
% nely: number of elements in y direction
% nelz: number of elements in z direction
% er: volume evolution rate
% volfrac: volume constraint
% E: Young's modulus, nu: Poisson’s ratio
% f: magnitude of load
% BESO3D(40,2,20,20,2,10,0.5,1,0.3,0.01,1.5) % example
[con,coord]=concoord3D(L,W,h,nelx,nely,nelz); % Element connectivity and Nodes' coordinate matrix
nel=size(con,1);nnd=size(coord,1); %number of nodes and elements
% Element elasticity matrix
D = (E/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0 ; NU 1-NU NU 0 0 0 ;
NU NU ...
1- NU 0 0 0 ;
0 0 0 (1-2*NU)/2 0 0 ; 0 0 0 0 (1- 2*NU)/2 0 ; 0 0 0 0 0 (1-2*NU)/2];
% initial design
x=ones(nel,1); % set all design variables as 1 (full solid design domain)
% Display and save the start mesh (full solid design domain)
surfcon=[];
for i=1:nel
surfcon(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5]);
con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
end
figure (1)
for i=1:size(surfcon,1)
for j=1:4
xcoord(j,i)=coord(surfcon(i,j),1);
ycoord(j,i)=coord(surfcon(i,j),2);
zcoord(j,i)=coord(surfcon(i,j),3);
end
end
h=patch(xcoord,ycoord,zcoord,'red');axis equal; axis([min(coord(:,1)) max(coord(:,1)) min(coord(:,2)) max(coord(:,2)) min(coord(:,3)) max(coord(:,3))]);axis off;
set(h,'LineWidth',1.5)
% set(h,'FaceColor',[.85 .85 .85])
set(h,'EdgeColor','red');
hgsave('0')
% Force and Boundary Conditions
U=sparse(3*nnd,1);
F = sparse(3*nnd,1);
fixeddofs=[1:3*(nely+1)*(nelz+1)]; % Fixed DoF
% fdof=3*(nelx*(nely+1)*(nelz+1) + (nely/2+1)*(nelz+1))-1; % DoF of the buttom of the free end
fdof=3*((nelx+1)*(nely+1)*(nelz+1)-(nely)*(nelz/2+1))-1;
F(fdof,1) = -1; %*nely/(nely+1); % apply force
alldofs = [1:3*nnd];
freedofs = setdiff(alldofs,fixeddofs);
% elements connected to each node
ncon=zeros(nnd,8);
for i=1:nnd
row=[];col=[];
[row,col] = find(con==i);
d=size(row,1);
ncon(i,1:d)=row;
end
%volume and stiffness of the elements
for el=1:nel
elcoord=coord(con(el,:)',:);
[~,Ae1(el,1)] = convhulln([elcoord(1,:);elcoord(2,:);elcoord(3,:);...
elcoord(4,:);elcoord(5,:);elcoord(6,:);elcoord(7,:);elcoord(8,:)]);
Ke1(:,:,el)=Brick_Stiffness(D,el,coord,con); % element stiffness at start mesh
end
A=sum(sum(Ae1)); % area of the whole structure
Ke=Ke1;Ae=Ae1; % set the current stifnees and volume as those of the initial design
if rmin>=0
[nbnodes,nbdis]=nodedis(coord,rmin);
end
volfracit=1; % defining initial volume fraction of the design
it=0;
while (((volfracit-volfrac)/volfrac)>=0.001)
    it=it+1;
if it>1; olddc=dc; end
volfracit=max(volfracit*(1-er),volfrac); % calculate volume of the next design
[K]=stiffness(nel,nnd,con,Ke,x); % Global stiffness matrix
U(freedofs,1)=K(freedofs,freedofs)\F(freedofs,1); % FEA for initial structure
for el = 1:nel
edof = [3*con(el,1)-2; 3*con(el,1)-1; 3*con(el,1)
3*con(el,2)-2; 3*con(el,2)-1; 3*con(el,2)
3*con(el,3)-2; 3*con(el,3)-1; 3*con(el,3)
3*con(el,4)-2; 3*con(el,4)-1; 3*con(el,4)
3*con(el,5)-2; 3*con(el,5)-1; 3*con(el,5)
3*con(el,6)-2; 3*con(el,6)-1; 3*con(el,6)
3*con(el,7)-2; 3*con(el,7)-1; 3*con(el,7)
3*con(el,8)-2; 3*con(el,8)-1; 3*con(el,8)];
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
if it>1; dc=(dc+olddc)/2; end % stabilization of evolutionary process
vf(it)=sum(sum(Ae))/A; % volume fraction of current design
% BESO Algorithm
[x]=ADDDEL(volfracit,A,dc,x,Ae);
sol=x==1;
vol(it)=sum(Ae1(sol,1));
%%%%%%%%%%%%%%%%%%%%%%%%% visualization of the topology
surfcon=[];xcoord=[];ycoord=[];zcoord=[];
x2coord=[];y2coord=[];z2coord=[];
surfcon2=[];
for i=1:nel
    if x(i)==1
        surfcon(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5]);
                con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
    else
        surfcon2(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5]);
                con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
    end
end
for i=1:size(surfcon,1)
for j=1:4
xcoord(j,i)=coord(surfcon(i,j),1);
ycoord(j,i)=coord(surfcon(i,j),2);
zcoord(j,i)=coord(surfcon(i,j),3);
end
end
for i=1:size(surfcon2,1)
for j=1:4
x2coord(j,i)=coord(surfcon2(i,j),1);
y2coord(j,i)=coord(surfcon2(i,j),2);
z2coord(j,i)=coord(surfcon2(i,j),3);
end
end
% figure (it)
% for i=1:size(surfcon,1)
% for j=1:4
% xcoord(j,i)=coord(surfcon(i,j),1);
% ycoord(j,i)=coord(surfcon(i,j),2);
% zcoord(j,i)=coord(surfcon(i,j),3);
% end
% end
% h2=patch(xcoord,ycoord,zcoord,'red');axis equal; axis([min(coord(:,1)) max(coord(:,1)) min(coord(:,2)) max(coord(:,2)) min(coord(:,3)) max(coord(:,3))]);axis off;
% h3=patch(x2coord,y2coord,z2coord,'blue');
% set(h2,'LineWidth',1.5)
% % set(h2,'FaceColor',[.85 .85 .85])
% set(h2,'EdgeColor','red');
% set(h3,'LineWidth',1.5)
% % set(h2,'FaceColor',[.85 .85 .85])
% set(h3,'EdgeColor','blue');
% % camlight
% % colormap(gray(100))
% % lighting gouraud
% hgsave(num2str(it))
% if it<200
% close(it);
% end
disp([' It.: ' sprintf('%4i',it) ' Obj.: ' sprintf('%10.4f',SE) 'Volfrac.: ' sprintf('%6.3f',vol(it)/A )]);
end
for w=1:10
it=it+1;
if it>1; olddc=dc; end
volfracit=max(volfracit*(1-er),volfrac); % calculate volume of the next design
[K]=stiffness(nel,nnd,con,Ke,x); % Global stiffness matrix
U(freedofs,1)=K(freedofs,freedofs)\F(freedofs,1); % FEA for initial structure
for el = 1:nel
edof = [3*con(el,1)-2; 3*con(el,1)-1; 3*con(el,1)
3*con(el,2)-2; 3*con(el,2)-1; 3*con(el,2)
3*con(el,3)-2; 3*con(el,3)-1; 3*con(el,3)
3*con(el,4)-2; 3*con(el,4)-1; 3*con(el,4)
3*con(el,5)-2; 3*con(el,5)-1; 3*con(el,5)
3*con(el,6)-2; 3*con(el,6)-1; 3*con(el,6)
3*con(el,7)-2; 3*con(el,7)-1; 3*con(el,7)
3*con(el,8)-2; 3*con(el,8)-1; 3*con(el,8)];
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
if it>1; dc=(dc+olddc)/2; end % stabilization of evolutionary process
vf(it)=sum(sum(Ae))/A; % volume fraction of current design
% BESO Algorithm
[x]=ADDDEL(volfracit,A,dc,x,Ae);
sol=x==1;
vol(it)=sum(Ae1(sol,1));
%%%%%%%%%%%%%%%%%%%%%%%%% visualization of the topology
surfcon=[];xcoord=[];ycoord=[];zcoord=[];
x2coord=[];y2coord=[];z2coord=[];
surfcon2=[];
m=0;
z=0;
for i=1:nel
    if x(i)==1
        surfcon(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5]);
                con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
    else
        surfcon2(end+1:end+6,1:4)=[con(i,1:4);con(i,5:8);con(i,[1 4 8 5]);
                con(i,[2 3 7 6]);con(i,[4 3 7 8]);con(i,[1 2 6 5])];
    end
end
for i=1:size(surfcon,1)
for j=1:4
xcoord(j,i)=coord(surfcon(i,j),1);
ycoord(j,i)=coord(surfcon(i,j),2);
zcoord(j,i)=coord(surfcon(i,j),3);
end
end
for i=1:size(surfcon2,1)
for j=1:4
x2coord(j,i)=coord(surfcon2(i,j),1);
y2coord(j,i)=coord(surfcon2(i,j),2);
z2coord(j,i)=coord(surfcon2(i,j),3);
end
end
figure (it)
% for i=1:size(surfcon,1)
% for j=1:4
% xcoord(j,i)=coord(surfcon(i,j),1);
% ycoord(j,i)=coord(surfcon(i,j),2);
% zcoord(j,i)=coord(surfcon(i,j),3);
% end
% end
h2=patch(xcoord,ycoord,zcoord,'red');axis equal; axis([min(coord(:,1)) max(coord(:,1)) min(coord(:,2)) max(coord(:,2)) min(coord(:,3)) max(coord(:,3))]);axis off;
% h3=patch(x2coord,y2coord,z2coord,'blue');
set(h2,'LineWidth',1.5)
% set(h2,'FaceColor',[.85 .85 .85])
set(h2,'EdgeColor','red');
% set(h3,'LineWidth',1.5)
% % set(h2,'FaceColor',[.85 .85 .85])
% set(h3,'EdgeColor','blue');
% camlight
% colormap(gray(100))
% lighting gouraud
hgsave(num2str(it))
% if it<200
% close(it);
% end

disp([' It.: ' sprintf('%4i',it) ' Obj.: ' sprintf('%10.4f',SE) 'Volfrac.: ' sprintf('%6.3f',vol(it)/A )]);
end
[K]=stiffness(nel,nnd,con,Ke,x); % Global stiffness matrix
U(freedofs,1)=K(freedofs,freedofs)\F(freedofs,1); % FEA for initial structure




ElementsElasticities=zeros(nel,1);
for i=1:nel
    if x(i)==1
       ElementsElasticities(i)=x(i);
    end
end
dlmwrite('Cantilever_40_10_10_50toisekato_E20_v02.txt', ElementsElasticities, 'delimiter', '\t');
end