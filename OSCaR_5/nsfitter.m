function [cplout] = nsfitter(cplin, const, kpin,F_ogin)
ns_neigh_dist=60;
%function to reduce the delta between north and south row ends
ns_delta=const.row_delta/2-0.05; %buffer the max offset by 1/10th
%find first and last tpzf in each row
B=~isnan(cplin.tpzf);
Indices=arrayfun(@(x) find(B(x,:), 1, 'last'), 1:size(cplin.tpzf,1)); %find the index of the final non-NaN per row
xf = arrayfun(@(x,y) cplin.tpx(y,x), Indices, 1:size(cplin.tpx,1)); %final tpx
yf = arrayfun(@(x,y) cplin.tpy(y,x), Indices, 1:size(cplin.tpy,1)); %final tpy
zf = arrayfun(@(x,y) cplin.tpzf(y,x), Indices, 1:size(cplin.tpzf,1)); %final tpz

xs=[cplin.tpx(:,1);xf'];
ys=[cplin.tpy(:,1);yf'];
zs=[cplin.tpzf(:,1);zf'];
rsi=[xs,ys];
rsf=[xs,ys];
[Idx,D]=rangesearch(rsi,rsf,ns_neigh_dist); %search for all row neighbors within neighbor_distance
maxLengthCell=max(cellfun('size',Idx,2));  %finding the longest vector in the cell array
for j=1:length(Idx)
    for k=cellfun('size',Idx(j),2)+1:maxLengthCell
        Idx{j}(k)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
    end
end
neighbors=cell2mat(Idx); %turning IDX into numerical values
neighbors(neighbors==0)=NaN; %turning 0s into NaN's so we can avoid working with them
[m,n]=size(neighbors);
nsnbr.tpx=xs;        nsnbr.tpy=ys;      nsnbr.tpz=zs;
nsnbr.n=zeros(m,1);  nsnbr.s=zeros(m,1);

for j=1:m
    for k=2:n
        if ~isnan(neighbors(j,k))
            if nsnbr.tpy(neighbors(j,k))<nsnbr.tpy(neighbors(j,1)) && abs(nsnbr.tpx(neighbors(j,k))-nsnbr.tpx(neighbors(j,1)))<5 %finding neigbors to south
                nsnbr.s(j)=neighbors(j,k);
            elseif nsnbr.tpy(neighbors(j,k))>nsnbr.tpy(neighbors(j,1)) && abs(nsnbr.tpx(neighbors(j,k))-nsnbr.tpx(neighbors(j,1)))<5 % finding neighbors to north
                nsnbr.n(j)=neighbors(j,k);
            end
        end
    end
end
nsnbr.nsd=zeros(length(nsnbr.tpx),1); %create zero array
for j=1:length(nsnbr.tpx) % writing this loop so north above south is always a positive number
    if ~nsnbr.n(j)==0
        if abs(nsnbr.tpz(j)-nsnbr.tpz(nsnbr.n(j)))>const.row_delta %if south minus north is greater than row delta (south above north)
            nsnbr.nsd(j,1)=nsnbr.tpz(j)-nsnbr.tpz(nsnbr.n(j)); %nsd = delta
        end
    elseif ~nsnbr.s(j)==0
        if abs(nsnbr.tpz(j)-nsnbr.tpz(nsnbr.s(j)))>const.row_delta %if north minus south is greater than row delta (north above south)
            nsnbr.nsd(j,1)=nsnbr.tpz(j)-nsnbr.tpz(nsnbr.s(j)); %nsd = delta
        end
    end
end
%nsnbr.nsd(nsnbr.nsd==0)=NaN; %turn zeros into NaN's
nsdel_orig=reshape(nsnbr.nsd,[],2); % full fledged deltas north and south (should be a matching *-1 for each)
nsdel_orig_n=nsdel_orig(:,1); nsdel_orig_s=nsdel_orig(:,2); %column 1 north, column 2 south
nsnbr.nsd(nsnbr.nsd>const.row_delta)=ns_delta; nsnbr.nsd(nsnbr.nsd<-const.row_delta)=-ns_delta; %creating target deltas
nsdel=reshape(nsnbr.nsd,[],2); %column 1 north, column 2 south %split back to initial and final
nsdel_n=nsdel(:,1); nsdel_s=nsdel(:,2);
nsdel_n_fin=(nsdel_orig_n)/2-nsdel_n;
nsdel_s_fin=(nsdel_orig_s)/2-nsdel_s;

calc_zi=zeros(1,length(zf)); calc_zf=zeros(1,length(zf));
nsdel_n_fin(nsdel_n_fin==0)=NaN; nsdel_s_fin(nsdel_s_fin==0)=NaN;

for i=1:length(zf)
    if ~isnan(nsdel_n_fin(i))
        calc_zi(i)=cplin.tpzf(i,1)-nsdel_n_fin(i);
    else
        calc_zi(i)=cplin.tpzf(i,1);
    end
    if ~isnan(nsdel_s_fin(i))
        calc_zf(i)=zf(i)-nsdel_s_fin(i);
    else
        calc_zf(i)=zf(i);
    end
end
pint=[cplin.tpx(:,1),cplin.tpy(:,1),calc_zi'];
pf=[xf',yf',calc_zf'];
uv=(pf-pint)./vecnorm(pf-pint,2,2); %find unit vector of all rows along each row of pf-pint (2 norm of dimension 2 [rows])
cplout.tpx=pint(:,1)-kpin.span.*uv(:,1);
cplout.tpy=pint(:,2)-kpin.span.*uv(:,2);
cplout.tpzf=pint(:,3)-kpin.span.*uv(:,3);
for i=1:length(Indices)
    cplout.hslp(i)=atand((calc_zf(i)-calc_zi(i))./abs(cplout.tpy((i),1)-cplout.tpy(i,Indices(i))));
end
cplout.hslp=cplout.hslp';
[~,n]=size(cplout.tpx);
cplout.slp=repmat(cplout.hslp,1,n);
cplout.bpz=F_ogin(cplout.tpx,cplout.tpy); %find z plumb to earth from top of cpl
cplout.bpy=cplout.tpy;
cplout.bpx=cplout.tpx;
cplout.mh=cplout.tpzf-cplout.bpz; %find length of all trusses
% % figure
% % scatter3(cplin.tpx(:,1),cplin.tpy(:,1),calc_zi,'ro','filled')
% % hold on
% % scatter3(cplin.tpx(:,1),cplin.tpy(:,1),cplin.tpzf(:,1),'rx')
% % scatter3(xf,yf,zf,'gx')
% % scatter3(xf,yf,calc_zf,'go','filled')
% % %view(2)
% % view(0,0)
% % nsnbr_table=struct2table(nsnbr);
% % ns_ugly=rmmissing(nsnbr_table);
% figure
% scatter3(reshape(cplout.tpx,[],1),reshape(cplout.tpy,[],1),reshape(cplout.tpzf,[],1));
% foo=1;
end