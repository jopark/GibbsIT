function emu=emufy(isotopomer,element)

l=log2(length(isotopomer));
nc=sum(de2bi(0:2^l-1),2);

if nargin==1 || isempty(element)
    m=l+1;
else
    m=min(length(element)+1,l+1);
    for i=1:l
        if min(i~=element)==1
            pat=2^i;
            nonelement=mod(1:length(nc),pat)';
            nonelement(nonelement==0)=pat;
            nonelement(nonelement<pat/2+1)=0;
            nonelement(nonelement>pat/2)=1;
            nc=nc-nonelement;
        end
    end
end
emu=zeros(1,m);
for i=1:m
    emu(i)=sum(isotopomer(nc==i-1));
end