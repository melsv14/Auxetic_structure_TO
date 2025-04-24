%%%%%%%%%% BESO add-remove algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=ADDDEL(volfrac,domain_vol,dc,x,Ae)
l1 = min(min(min(dc))); l2 = max(max(max(dc)));
while ((l2-l1)/l2 > 1.0e-5)
th = (l1+l2)/2.0;
x = max(0.0001,sign(dc-th));
z=0;
sizeEl=size(dc,1);
for i=1:size(dc,1)
    if x(i)==1
        z=z+1;
    end
end
if z*Ae-volfrac*(domain_vol) > 0
l1 = th;
else
l2 = th;
end
end