function [mname] = sum_of_sums(const, kpv, pd_out)
for i=1:length(kpv)
    sos.OBJECTID=pd_out.(kpv{i}).bsr;
    sos.BPX=pd_out.(kpv{i}).bpx;
    sos.BPY=pd_out.(kpv{i}).bpy;
    for j=1:length(pd_out.(kpv{i}).pt_desc)
        parts=strsplit(string(pd_out.(kpv{i}).pt_desc(j)),'_');
        p1{j}=parts{1};
        p2{j}=parts{2};
    end
    bpf=find(strcmpi(string(p2),"std"));
    gpf=find(strcmpi(string(p2),"mtr"));
    sos.BPGP=cell(length(pd_out.(kpv{i}).pt_desc),1);
    sos.BPGP(bpf)={'BP'};
    sos.BPGP(gpf)={'GP'};
    sos.EXTINT=cell(length(pd_out.(kpv{i}).pt_desc),1);
    ipf=find(strcmpi(string(p1),"int"));
    epf=find(strcmpi(string(p1),"ext"));
    sos.EXTINT(ipf)={"INT"};
    sos.EXTINT(epf)={"EXT"};
    sos.TTH=pd_out.(kpv{i}).pilebin;
    sosout=struct2table(sos);
    if const.writefiles==1
        for i=1:length(kpv)
            if const.simple_out==0
                mkdir([const.outpath '/' kpv{i} '_output']);
                writetable(sosout,[const.outpath '/' kpv{i} '_output/' kpv{i} '_SumofSums_input.xlsx']);
                mname=[const.outpath '/' kpv{i} '_output/' kpv{i} '_SumofSums_input.xlsx'];
            else
                mkdir([const.outpath '/' kpv{i} '_SumofSums']);
                writetable(sosout,[const.outpath '/' kpv{i} '_SumofSums/' kpv{i} '_SumofSums_input.xlsx']);
                mname=[const.outpath '/' kpv{i} '_SumofSums/' kpv{i} '_SumofSums_input.xlsx'];
            end
        end
    end
end
end
