function[ke,globalGP]=stiffnessmat(D,en,coord,con)
%%%%%%%% Element Stiffness Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GP=[-1/sqrt(3),-1/sqrt(3);1/sqrt(3),-1/sqrt(3);1/sqrt(3),1/sqrt(3);-1/sqrt(3),1/sqrt(3)];
i=con(en,1);j=con(en,2);k=con(en,3);l=con(en,4);
x01=coord(i,1);y01=coord(i,2);
x02=coord(j,1);y02=coord(j,2);
x03=coord(k,1);y03=coord(k,2);
x04=coord(l,1);y04=coord(l,2);
%for the visualization of gp
coordinatesOfElement=[x01 y01;x02 y02;x03 y03;x04 y04];
globalGP=zeros(4,2);
ke=zeros(8,8);
for i=1:4
r=GP(i,1);s=GP(i,2);
N1=1/4*(1-r)*(1-s);
N2=1/4*(1+r)*(1-s);
N3=1/4*(1+r)*(1+s);
N4=1/4*(1-r)*(1+s);
N=[N1;N2;N3;N4];
%for the visualization of gp
globalGP(i,:)=transpose(N)*coordinatesOfElement;
dN1r=-1/4*(1-s);dN1s=-1/4*(1-r);
dN2r=1/4*(1-s);dN2s=-1/4*(1+r);
dN3r=1/4*(1+s);dN3s=1/4*(1+r);
dN4r=-1/4*(1+s);dN4s=1/4*(1-r);
j11=x01*dN1r+x02*dN2r+x03*dN3r+x04*dN4r;
j12=y01*dN1r+y02*dN2r+y03*dN3r+y04*dN4r;
j21=x01*dN1s+x02*dN2s+x03*dN3s+x04*dN4s;
j22=y01*dN1s+y02*dN2s+y03*dN3s+y04*dN4s;
J=[j11 j12;j21 j22];
detj=det(J);
dNxy=inv(J)*[dN1r dN2r dN3r dN4r; dN1s dN2s dN3s dN4s];
dN1x=dNxy(1,1);
dN2x=dNxy(1,2);
dN3x=dNxy(1,3);
dN4x=dNxy(1,4);
dN1y=dNxy(2,1);
dN2y=dNxy(2,2);
dN3y=dNxy(2,3);
dN4y=dNxy(2,4); %the ; were missed in BLO below...!!!!
BL0=[dN1x 0 dN2x 0 dN3x 0 dN4x 0;
    0 dN1y 0 dN2y 0 dN3y 0 dN4y;
    dN1y dN1x dN2y dN2x dN3y dN3x dN4y dN4x];
ke=ke+BL0'*D*BL0*detj;
end