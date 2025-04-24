%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates element stifness matrix for 8-node Hex elements.
% by Meisam Abdi, university of Nottinham
function [ke] = Brick_Stiffness(D,en,coord,con)
g=1/sqrt(3);
GP=[-g -g -g;g -g -g;g g -g;-g g -g;-g -g g;g -g g;g g g;-g g g];
%integration points
i=con(en,1);
j=con(en,2);
k=con(en,3);
l=con(en,4);
m=con(en,5);
n=con(en,6);
o=con(en,7);
p=con(en,8);
x1=coord(i,1);y1=coord(i,2);z1=coord(i,3);
x2=coord(j,1);y2=coord(j,2);z2=coord(j,3);
x3=coord(k,1);y3=coord(k,2);z3=coord(k,3);
x4=coord(l,1);y4=coord(l,2);z4=coord(l,3);
x5=coord(m,1);y5=coord(m,2);z5=coord(m,3);
x6=coord(n,1);y6=coord(n,2);z6=coord(n,3);
x7=coord(o,1);y7=coord(o,2);z7=coord(o,3);
x8=coord(p,1);y8=coord(p,2);z8=coord(p,3);
% shape functions
% N1 = 1/8*(1-s)*(1-t)*(1-u); N2 = 1/8*(1+s)*(1-t)*(1-u);
% N3 = 1/8*(1+s)*(1+t)*(1-u); N4 = 1/8*(1-s)*(1+t)*(1-u);
% N5 = 1/8*(1-s)*(1-t)*(1+u); N6 = 1/8*(1+s)*(1-t)*(1+u);
% N7 = 1/8*(1+s)*(1+t)*(1+u); N8 = 1/8*(1-s)*(1+t)*(1+u);
% interpolation of shape functions
% x = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 +N8*x8;
% y = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 +N8*y8;
% z = N1*z1 + N2*z2 + N3*z3 + N4*z4 + N5*z5 + N6*z6 + N7*z7 +N8*z8;
ke=zeros(24,24);
for ii=1:8
    s=GP(ii,1);t=GP(ii,2);u=GP(ii,3);
xs = (x2*(t - 1)*(u - 1))/8 - (x1*(t - 1)*(u - 1))/8 - (x3*(t +1)*(u - 1))/8 + (x4*(t + 1)*(u - 1))/8 + (x5*(t - 1)*(u + 1))/8-(x6*(t - 1)*(u + 1))/8 + (x7*(t + 1)*(u + 1))/8 - (x8*(t + 1)*(u+1))/8; 
% diff(x,t);
xt = x2*(s/8 + 1/8)*(u - 1) - x1*(s/8 - 1/8)*(u - 1) - x3*(s/8 +1/8)*(u - 1)+x4*(s/8 - 1/8)*(u - 1) + x5*(s/8 - 1/8)*(u + 1)-x6*(s/8 + 1/8)*(u + 1) + x7*(s/8 + 1/8)*(u + 1)-x8*(s/8 -1/8)*(u + 1); 
% diff(x,t);
xu = x2*(s/8 + 1/8)*(t - 1) - x1*(s/8 - 1/8)*(t - 1) - x3*(s/8+1/8)*(t + 1) + x4*(s/8 - 1/8)*(t + 1) + x5*(s/8 - 1/8)*(t - 1) -x6*(s/8 + 1/8)*(t - 1) + x7*(s/8 + 1/8)*(t + 1) - x8*(s/8 -1/8)*(t + 1);
% diff(x,u);
ys = (y2*(t - 1)*(u - 1))/8 - (y1*(t - 1)*(u - 1))/8 - (y3*(t+1)*(u - 1))/8 + (y4*(t + 1)*(u - 1))/8 + (y5*(t - 1)*(u + 1))/8-(y6*(t - 1)*(u + 1))/8 + (y7*(t + 1)*(u + 1))/8 - (y8*(t + 1)*(u+1))/8; 
% diff(y,s);
yt = y2*(s/8 + 1/8)*(u - 1) - y1*(s/8 - 1/8)*(u - 1) - y3*(s/8+1/8)*(u - 1) + y4*(s/8 - 1/8)*(u - 1) + y5*(s/8 - 1/8)*(u + 1)-y6*(s/8 + 1/8)*(u + 1) + y7*(s/8 + 1/8)*(u + 1) - y8*(s/8-1/8)*(u + 1); 
% diff(y,t);
yu = y2*(s/8+1/8)*(t - 1) - y1*(s/8 - 1/8)*(t - 1) - y3*(s/8+1/8)*(t + 1) + y4*(s/8 - 1/8)*(t + 1) + y5*(s/8 - 1/8)*(t - 1) -y6*(s/8 + 1/8)*(t - 1) + y7*(s/8 + 1/8)*(t + 1) - y8*(s/8-1/8)*(t + 1); 
% diff(y,u);
zs = (z2*(t - 1)*(u - 1))/8 - (z1*(t - 1)*(u - 1))/8 - (z3*(t+1)*(u - 1))/8 + (z4*(t + 1)*(u - 1))/8 + (z5*(t - 1)*(u + 1))/8-(z6*(t - 1)*(u + 1))/8 + (z7*(t + 1)*(u + 1))/8 - (z8*(t + 1)*(u+1))/8; 
% diff(z,s);
zt=z2*(s/8 + 1/8)*(u - 1) - z1*(s/8 - 1/8)*(u - 1) - z3*(s/8+1/8)*(u - 1) + z4*(s/8 - 1/8)*(u - 1) + z5*(s/8 - 1/8)*(u + 1)-z6*(s/8 + 1/8)*(u + 1) + z7*(s/8 + 1/8)*(u + 1) - z8*(s/8-1/8)*(u + 1); 
% diff(z,t);
zu=z2*(s/8 + 1/8)*(t - 1) - z1*(s/8 - 1/8)*(t - 1) - z3*(s/8+1/8)*(t + 1) + z4*(s/8 - 1/8)*(t + 1) + z5*(s/8 - 1/8)*(t - 1)-z6*(s/8 + 1/8)*(t - 1) + z7*(s/8 + 1/8)*(t + 1) - z8*(s/8-1/8)*(t + 1); 
% diff(z,u);
J = [xs ys zs; xt yt zt; xu yu zu];
detJ = xs*(yt*zu - zt*yu) - ys*(xt*zu - zt*xu) + zs*(xt*yu -yt*xu);
N1s=-((t - 1)*(u - 1))/8; % diff(N1,s);
N2s=((t - 1)*(u - 1))/8; % diff(N2,s);
N3s=-((t + 1)*(u - 1))/8; % diff(N3,s);
N4s=((t + 1)*(u - 1))/8; % diff(N4,s);
N5s=((t - 1)*(u + 1))/8; % diff(N5,s);
N6s=-((t - 1)*(u + 1))/8; % diff(N6,s);
N7s=((t + 1)*(u + 1))/8; % diff(N7,s);
N8s=-((t + 1)*(u + 1))/8; % diff(N8,s);
N1t = -(s/8 - 1/8)*(u - 1); % diff(N1,t);
N2t = (s/8 + 1/8)*(u - 1); % diff(N2,t);
N3t = -(s/8 + 1/8)*(u - 1); % diff(N3,t);
N4t = (s/8 - 1/8)*(u - 1); % diff(N4,t);
N5t = (s/8 - 1/8)*(u + 1); % diff(N5,t);
N6t = -(s/8 + 1/8)*(u + 1); % diff(N6,t);
N7t = (s/8 + 1/8)*(u + 1); % diff(N7,t);
N8t = -(s/8 - 1/8)*(u + 1); % diff(N8,t);
N1u = -(s/8 - 1/8)*(t - 1); % diff(N1,u);
N2u = (s/8 + 1/8)*(t - 1); % diff(N2,u);
N3u = -(s/8 + 1/8)*(t + 1); % diff(N3,u);
N4u = (s/8 - 1/8)*(t + 1); % diff(N4,u);
N5u = (s/8 - 1/8)*(t - 1); % diff(N5,u);
N6u = -(s/8 + 1/8)*(t - 1); %diff(N6,u);
N7u = (s/8 + 1/8)*(t + 1); % diff(N7,u);
N8u = -(s/8 - 1/8)*(t + 1); % diff(N8,u);
Nxyz=J\[N1s N2s N3s N4s N5s N6s N7s N8s;
N1t N2t N3t N4t N5t N6t N7t N8t;
N1u N2u N3u N4u N5u N6u N7u N8u];
N1x=Nxyz(1,1); N2x=Nxyz(1,2); N3x=Nxyz(1,3); N4x=Nxyz(1,4);
N5x = Nxyz(1,5); N6x = Nxyz(1,6); N7x = Nxyz(1,7); N8x =Nxyz(1,8);
N1y = Nxyz(2,1); N2y = Nxyz(2,2); N3y = Nxyz(2,3); N4y =Nxyz(2,4);
N5y = Nxyz(2,5); N6y = Nxyz(2,6); N7y = Nxyz(2,7); N8y =Nxyz(2,8);
N1z = Nxyz(3,1); N2z = Nxyz(3,2); N3z = Nxyz(3,3); N4z =Nxyz(3,4);
N5z = Nxyz(3,5); N6z = Nxyz(3,6); N7z = Nxyz(3,7); N8z =Nxyz(3,8);
% Linear strain-displacement transformation matrix
B=[N1x 0 0 N2x 0 0 N3x 0 0 N4x 0 0 N5x 0 0 N6x 0 0 N7x 0 0 N8x 0 0;
0 N1y 0 0 N2y 0 0 N3y 0 0 N4y 0 0 N5y 0 0 N6y 0 0 N7y 0 0 N8y 0;
0 0 N1z 0 0 N2z 0 0 N3z 0 0 N4z 0 0 N5z 0 0 N6z 0 0 N7z 0 0 N8z;
N1y N1x 0 N2y N2x 0 N3y N3x 0 N4y N4x 0 N5y N5x 0 N6y N6x 0 N7y N7x 0 N8y N8x 0;
0 N1z N1y 0 N2z N2y 0 N3z N3y 0 N4z N4y 0 N5z N5y 0 N6z N6y 0 N7z N7y 0 N8z N8y;
N1z 0 N1x N2z 0 N2x N3z 0 N3x N4z 0 N4x N5z 0 N5x N6z 0 N6x N7z 0 N7x N8z 0 N8x];
BD=transpose(B)*D*B*detJ;
ke=ke+BD;
end
