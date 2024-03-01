function [drow] = a_8_rowfits_6(const, surface, drow)
%SLOPE SOLUTION - uncorrected raw surface heights
drow.npz=surface.F_og(drow.npx,drow.npy); %find npz for initial fit
drow.spz=surface.F_og(drow.npx,drow.spy); %find spz for initial fit
for i=1:height(drow)
    xs=drow.npx(i);
    ys=drow.spy(i);
    yn=drow.npy(i);
    zs=drow.spz(i);
    zn=drow.npz(i);
    calchyp=hypot((yn-ys),(zn-zs));
    yd=cosd(atand((zn-zs)/(yn-ys)))*(calchyp-(yn-ys));
    ys=ys-yd;
    zs=surface.F_og(xs,ys);
    drow.spyc(i)=ys;
    drow.spzc(i)=zs;
end
drow.rlengthc=drow.npy-drow.spyc;
%ROW BEST FIT - based on shortened row above
for i=1:height(drow)
    drowh.yinter{i,1}=drow.npy(i):-const.bin_y:drow.spyc(i); %create query/polyfit points along each of the rows
    drowh.xinter{i,1}=drow.npx(i).*ones(1,length(drowh.yinter{i})); %assign x's for z determination
    drowh.zinter{i,1}=surface.F_og(drowh.xinter{i},drowh.yinter{i}); %determine z's at each of the y locations
    drowh.pf_c{i,1}=polyfit(drowh.yinter{i},drowh.zinter{i},1); %determine polyfit coefficients
    drowh.pf_f{i,1}=polyval(drowh.pf_c{i},drowh.yinter{i}); %line of linear best fit of the row
    drow.lbf_npz(i)=drowh.pf_f{i}(1); %linear best fit z of northern point
    drow.lbf_spz(i)=drowh.pf_f{i}(end); %linear best fit z of southern point
end

%E/W BEST FIT
drowh.ns_comb_x=[drow.npx, drow.npx]; %this assumes straight n/s rows
drowh.ns_comb_y=[drow.npy, drow.spyc];
drowh.ns_comb_z=[drow.npz, drow.spzc]; %this assumes straight n/s rows
allneighs=[drow.row,drow.nw_id,drow.ne_id]; %all rows within neighbor distance, including target row

for i=1:length(allneighs)
    hld=allneighs(i,:);
    anr{i,:}=hld(~isnan(hld)); %reduce allneighs to only non-nan rows (all neighs reduced = anr)
end

for i=1:numel(anr)
    pfn1=polyfit(drow.npx(anr{i}),drow.npz(anr{i}),2); %this is all set up assuming rows are perfect n/s
    pfs1=polyfit(drow.npx(anr{i}),drow.spzc(anr{i}),2);
    pfnz1{i}=polyval(pfn1,drow.npx(anr{i}));
    pfsz1{i}=polyval(pfs1,drow.npx(anr{i}));
    drow.nbf_npz(i)=pfnz1{i}(1);
    drow.nbf_spz(i)=pfsz1{i}(1);
end

drow=removevars(drow,{'nw_id','ne_id'});


%check out some rows
% tr=25;
% figure
% scatter(drow.npx(anr{tr}),drow.npz(anr{tr}),'red');
% hold on
% scatter(drow.npx(anr{tr}),drow.nbf_npz(anr{tr}),'green');
end