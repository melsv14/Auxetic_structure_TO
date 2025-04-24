%%%%%%%%%%%%% filteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b]=nodedis(coord,rmin)
for i=1:size(coord,1)
coordi=zeros(size(coord));
coordi(:,1)=coord(i,1);coordi(:,2)=coord(i,2);
ndis=coordi-coord;
for j=1:size(coord,1)
r(i,j)=sqrt(sum(ndis(j,1)^2+ndis(j,2)^2));
end
a{i}=find(r(i,:)<rmin);
b{i}=r(i,a{i});
end
end