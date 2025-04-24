function [ndc,dc]=filterBESO(a,b,ndc,dc,rmin,con)
oldndc=ndc;
for i=1:size(ndc,1)
wi=rmin-b{i};
alphai=oldndc(a{i},1)';
ndc(i,1)=sum(wi.*alphai)/sum(wi);
end
for i=1:size(dc,1)
dc(i,1)=mean(ndc(con(i,:)',1));
end
end