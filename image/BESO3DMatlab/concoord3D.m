function [connectivity,coordinate]=concoord3D(L,W,h,nelx,nely,nelz)
connectivity=zeros(nelx*nely*nelz,8);
coordinate=zeros((nelx+1)*(nely+1)*(nelz+1),3);
for ii=1:nelx*nely*nelz
    elx=fix((ii-1)/(nely*nelz))+1;
    ely=fix(((ii-(elx-1)*nely*nelz)-1)/nelz)+1;
    elz=mod(ii,nelz);
    if elz==0
        elz=nelz;
    end
    connectivity(ii,1)=(nely+1)*(nelz+1)*(elx-1)+(ely-1)*(nelz+1)+elz;
    connectivity(ii,2)=connectivity(ii,1)+(nely+1)*(nelz+1);
    connectivity(ii,3)=connectivity(ii,2)+(nelz+1);
    connectivity(ii,4)=connectivity(ii,1)+(nelz+1);
    connectivity(ii,5)=connectivity(ii,1)+1;
    connectivity(ii,6)=connectivity(ii,2)+1;
    connectivity(ii,7)=connectivity(ii,3)+1;
    connectivity(ii,8)=connectivity(ii,4)+1;
end
dx=L/nelx;
% dy=W/nely;
% dz=h/nelz;
dy=h/nely;
dz=W/nelz;
for ii=1:(nelx+1)*(nely+1)*(nelz+1)
       nx=fix((ii-1)/((nely+1)*(nelz+1)))+1;
       ny=fix(((ii-(nx-1)*(nely+1)*(nelz+1))-1)/(nelz+1))+1;
       nz=mod(ii,(nelz+1));
       if nz==0
           nz=nelz+1;
       end
coordinate(ii,:)=[(nx-1)*dx (ny-1)*dy (nz-1)*dz];
end
end