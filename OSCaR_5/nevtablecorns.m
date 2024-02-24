function nevtablecorns(cpl, kpv, pmd, const)
pan_len=const.xpdim/25.4/12;
for i=1:length(kpv)
    daly.sect=pmd.(kpv{i}).sect;
    daly.row=pmd.(kpv{i}).row;
    daly.tpx=reshape(cpl.(kpv{i}).tpx,[],1);
    daly.tpy=reshape(cpl.(kpv{i}).tpy,[],1);
    daly.tpz=reshape(cpl.(kpv{i}).tpzf,[],1);
    daly.upslope=reshape(cpl.(kpv{i}).slopeu',[],1);
    daly.downslope=reshape(cpl.(kpv{i}).sloped',[],1);
    daly.upslope(daly.upslope==0)=NaN;
    daly.downslope(daly.downslope==0)=NaN;
    daly.westx=daly.tpx-pan_len/2;
    daly.eastx=daly.tpx+pan_len/2;
    daly.westy=daly.tpy;
    daly.easty=daly.tpy;
    daly.westz=daly.tpz+0.25;
    daly.eastz=daly.tpz+0.25;

    dalyout=struct2table(daly);
    dalyout=rmmissing(dalyout);
    dalyout.(kpv{i})=dalyout;
end

for i=1:length(kpv)
    if const.simple_out==0
        mkdir([const.outpath '/' kpv{i} '_DalyPoints']);
        writetable(dalyout.(kpv{i}),[const.outpath '/' kpv{i} '_DalyPoints/' kpv{i} '_DalyPts.csv']);
        %writetable(ld_out.(kpv{i}),[const.outpath '/legdata/' kpv{i} '/' 'legData.csv']);
    else
        mkdir([const.outpath '/' kpv{i} '_DalyPoints']);
        writetable(dalyout.(kpv{i}),[const.outpath '/' kpv{i} '_DalyPoints/' kpv{i} '_DalyPts.csv']);
    end
end

end