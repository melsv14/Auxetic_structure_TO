%%%%%%%%%% global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=stiffness(nel,nnd,con,Ke,x)
K = sparse(3*nnd, 3*nnd);
for el = 1:nel
edof = [3*con(el,1)-2; 3*con(el,1)-1; 3*con(el,1)
3*con(el,2)-2; 3*con(el,2)-1; 3*con(el,2)
3*con(el,3)-2; 3*con(el,3)-1; 3*con(el,3)
3*con(el,4)-2; 3*con(el,4)-1; 3*con(el,4)
3*con(el,5)-2; 3*con(el,5)-1; 3*con(el,5)
3*con(el,6)-2; 3*con(el,6)-1; 3*con(el,6)
3*con(el,7)-2; 3*con(el,7)-1; 3*con(el,7)
3*con(el,8)-2; 3*con(el,8)-1; 3*con(el,8)];
K(edof,edof) = K(edof,edof) + Ke(:,:,el)*x(el);
end
%%%%%%%%%%%%%%
end