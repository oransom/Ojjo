function [ns_ugly] = northsouth(kpv, cpl, const)
ns_delta=const.row_delta;
neighbor_distance=40;
ns_all=cell2table(cell(0,4), 'VariableNames', {'tpx','tpy','tpz','pt'});
for i=1:length(kpv)
    ns.(kpv{i}).tpx=reshape(cpl.(kpv{i}).tpx,[],1);
    ns.(kpv{i}).tpy=reshape(cpl.(kpv{i}).tpy,[],1);
    ns.(kpv{i}).tpz=reshape(cpl.(kpv{i}).tpzf,[],1);
    ns.(kpv{i}).pt=reshape(cpl.(kpv{i}).piletype,[],1);
    ns.(kpv{i}).pt(ns.(kpv{i}).pt==0)=NaN;
    nst.(kpv{i})=struct2table(ns.(kpv{i}));
    nst.(kpv{i})=rmmissing(nst.(kpv{i}));
    ns_all=[ns_all; nst.(kpv{i})];
end


rsx=[ns_all.tpx, ns_all.tpy];
rsy=[ns_all.tpx, ns_all.tpy];
[Idx,D]=rangesearch(rsx,rsy,neighbor_distance); %search for all row neighbors within neighbor_distance
maxLengthCell=max(cellfun('size',Idx,2));  %finding the longest vector in the cell array
for i=1:length(Idx)
   for j=cellfun('size',Idx(i),2)+1:maxLengthCell
       Idx{i}(j)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
   end
end
neighbors=cell2mat(Idx); %turning IDX into numerical values
neighbors(neighbors==0)=NaN; %turning 0s into NaN's so we can avoid working with them

[m,n]=size(neighbors);
nsnbr.tpx=ns_all.tpx;        nsnbr.tpy=ns_all.tpy;      nsnbr.tpz=ns_all.tpz;     nsnbr.pt=ns_all.pt;
nsnbr.n=zeros(m,1);          nsnbr.s=zeros(m,1);   
%nsnbr.bsre=cell(m,1);        nsnbr.movecr=zeros(m,1);   nsnbr.canmove=zeros(m,1);       %allocate matrices for row neighbor information

for i=1:m
    for j=2:n
        if ~isnan(neighbors(i,j))
            if nsnbr.tpy(neighbors(i,j))<nsnbr.tpy(neighbors(i,1)) && abs(nsnbr.tpx(neighbors(i,j))-nsnbr.tpx(neighbors(i,1)))<3 %finding neigbors to west
                nsnbr.s(i)=neighbors(i,j);
            elseif nsnbr.tpy(neighbors(i,j))>nsnbr.tpy(neighbors(i,1)) && abs(nsnbr.tpx(neighbors(i,j))-nsnbr.tpx(neighbors(i,1)))<3 % finding neighbors to east
                nsnbr.n(i)=neighbors(i,j);
            end
        end
    end
end
nsnbr.nsd=zeros(length(nsnbr.tpx),1); %create zero array
for i=1:length(nsnbr.tpx)
    if ~nsnbr.n(i)==0
        if abs(nsnbr.tpz(i)-nsnbr.tpz(nsnbr.n(i)))>ns_delta
            nsnbr.nsd(i,1)=nsnbr.tpz(i)-nsnbr.tpz(nsnbr.n(i));
        end
    elseif ~nsnbr.s(i)==0
        if abs(nsnbr.tpz(i)-nsnbr.tpz(nsnbr.s(i)))>ns_delta
            nsnbr.nsd(i,1)=nsnbr.tpz(i)-nsnbr.tpz(nsnbr.s(i));
        end
    end
end
nsnbr.nsd(nsnbr.nsd==0)=NaN;
nsnbr_table=struct2table(nsnbr);
ns_ugly=rmmissing(nsnbr_table);
end
