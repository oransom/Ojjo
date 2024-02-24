function [ns_lim_out] = nslimiter(cplin, const, kpin, F_ogin)
%function to limit north and south row slope to allowable levels
%find first and last tpzf in each row
B=~isnan(cplin.tpzf);
Indices=arrayfun(@(x) find(B(x,:), 1, 'last'), 1:size(cplin.tpzf,1)); %find the index of the final non-NaN per row
xf = arrayfun(@(x,y) cplin.tpx(y,x), Indices, 1:size(cplin.tpx,1))'; %final tpx
yf = arrayfun(@(x,y) cplin.tpy(y,x), Indices, 1:size(cplin.tpy,1))'; %final tpy
zf = arrayfun(@(x,y) cplin.tpzf(y,x), Indices, 1:size(cplin.tpzf,1))'; %final tpz

xi=[cplin.tpx(:,1)];
yi=[cplin.tpy(:,1)];
zi=[cplin.tpzf(:,1)];

calc_slp=-atand((zi-zf)./(yi-yf));
idn=find(calc_slp>0 & calc_slp>const.nslopelimit);
ids=find(calc_slp<0 & calc_slp<-const.sslopelimit);

ns_diff=calc_slp(idn)-const.nslopelimit;
ss_diff=calc_slp(ids)+const.sslopelimit;
ns_ydiff=yi(idn)-yf(idn);
ss_ydiff=yi(ids)-yf(ids);

zn_diff_calc=tand(ns_diff).*ns_ydiff;
zs_diff_calc=tand(ss_diff).*ss_ydiff;

for i=1:length(zn_diff_calc)
    if zi(idn(i))>zf(idn(i))
        zi(idn(i))=zi(idn(i))-zn_diff_calc(i)/2;
        zf(idn(i))=zf(idn(i))+zn_diff_calc(i)/2;
    else
        zi(idn(i))=zi(idn(i))+zn_diff_calc(i)/2;
        zf(idn(i))=zf(idn(i))-zn_diff_calc(i)/2;
    end
end
for i=1:length(zs_diff_calc)
    if zi(ids(i))>zf(ids(i))
        zi(ids(i))=zi(ids(i))+zs_diff_calc(i)/2;
        zf(ids(i))=zf(ids(i))-zs_diff_calc(i)/2;
    else
        zi(ids(i))=zi(ids(i))-zs_diff_calc(i)/2;
        zf(ids(i))=zf(ids(i))+zs_diff_calc(i)/2;
    end
end

calc_slp_chk=-atand((zi-zf)./(yi-yf));

pint=[xi, yi, zi];
pf=[xf, yf, zf];
uv=(pf-pint)./vecnorm(pf-pint,2,2); %find unit vector of all rows along each row of pf-pint (2 norm of dimension 2 [rows])
ns_lim_out.tpx=pint(:,1)-kpin.span.*uv(:,1);
ns_lim_out.tpy=pint(:,2)-kpin.span.*uv(:,2);
ns_lim_out.tpzf=pint(:,3)-kpin.span.*uv(:,3);
[~,n]=size(ns_lim_out.tpx);
ns_lim_out.slp=repmat(calc_slp_chk,1,n);
ns_lim_out.bpz=F_ogin(ns_lim_out.tpx,ns_lim_out.tpy); %find z plumb to earth from top of cpl
ns_lim_out.bpy=ns_lim_out.tpy;
ns_lim_out.bpx=ns_lim_out.tpx;
ns_lim_out.mh=ns_lim_out.tpzf-ns_lim_out.bpz; %find length of all trusses
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
% scatter3(reshape(ns_lim_out.tpx,[],1),reshape(ns_lim_out.tpy,[],1),reshape(ns_lim_out.tpzf,[],1));
%foo=1;
end