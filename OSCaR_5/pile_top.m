function [top_out] = pile_top(kpv, pmd, po, ppp_out)

for i=1:length(kpv)
    %create pinning table
    %pin.(kpv{i}).bsr=pmd.(kpv{i}).bsr';
    top.(kpv{i}).p=pmd.(kpv{i}).bsrt';
    top.(kpv{i}).n=po.(kpv{i}).tpy;
    top.(kpv{i}).e=po.(kpv{i}).tpx;
    top.(kpv{i}).z=po.(kpv{i}).tpz;
    top.(kpv{i}).d=pmd.(kpv{i}).tt;
    top_out.(kpv{i})=struct2table(top.(kpv{i}));
    top_out.(kpv{i})=rmmissing(top_out.(kpv{i}));
    for j=1:length(top_out.(kpv{i}).z)
        top_out.(kpv{i}).z(j)=ppp_out.(kpv{i}).top_pile_elev_ft(j);
    end
end
