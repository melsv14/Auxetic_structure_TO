%%%%%%%%%% global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=stiffness(nel,nnd,con,Ke,x)
% K = zeros(2*nnd, 2*nnd);
K = sparse(2*nnd, 2*nnd);
for el = 1:nel
% Boole=zeros(8,2*nnd);   
edof = [2*con(el,1)-1; 2*con(el,1); 2*con(el,2)-1;
2*con(el,2); 2*con(el,3)-1; 2*con(el,3); 2*con(el,4)-1;
2*con(el,4)];
% for i=1:8
%     Boole(i,edof(i))=1;
% end
%     K=K+transpose(Boole)*Ke(:,:,el)*Boole;
K(edof,edof) = K(edof,edof) + Ke(:,:,el)*x(el);
end
end
