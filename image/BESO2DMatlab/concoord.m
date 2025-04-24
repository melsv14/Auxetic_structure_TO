function [connectivity,coord]=concoord(nelx,nely,L,h)
connectivity=zeros(nelx*nely,4);
for ii=1:nelx*nely
  rw=mod(ii,nely);
  cl=fix((ii-1)/nely)+1;
  connectivity(ii,1)=cl-1+ii;
  connectivity(ii,2)=connectivity(ii,1)+nely+1;
  connectivity(ii,3)=connectivity(ii,2)+1;
  %the code has not values for the 4th column!!!
  connectivity(ii,4)=connectivity(ii,1)+1;
end
coord=zeros((nelx+1)*(nely+1),2);
for ii=1:(nelx+1)*(nely+1)
coord(ii,1)=fix((ii-1)/(nely+1))*L/nelx;
coord(ii,2)=mod((ii-1),(nely+1))*h/nely;
end
end