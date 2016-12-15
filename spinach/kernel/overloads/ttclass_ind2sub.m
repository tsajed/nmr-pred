function ivec=ttclass_ind2sub(siz,ind)
d=numel(siz);
ind=ind-1;
ivec=zeros(d,1);
for k=d:-1:1
    ivec(k)=mod(ind,siz(k));
    ind=floor(ind/siz(k));
end
ivec=ivec+1;
end