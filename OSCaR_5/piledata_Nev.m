function [pdout, ot, pileTable] = piledata_Nev(kpv, pmd, po, to, oos_data, const)
%geo=readtable(string(const.geo_inputs));
lf=const.leg_factory{1};
for i=1:length(kpv)
    pd.(kpv{i}).bsr=pmd.(kpv{i}).bsr'; %just dumping the block section row identifier in here to make it work for now
    pd.(kpv{i}).block=pmd.(kpv{i}).block';
    pd.(kpv{i}).section=pmd.(kpv{i}).sect;
    pd.(kpv{i}).row_number=pmd.(kpv{i}).row;
    pd.(kpv{i}).truss_letter=pmd.(kpv{i}).tl;
    pd.(kpv{i}).tpx=po.(kpv{i}).tpx;
    pd.(kpv{i}).tpy=po.(kpv{i}).tpy;
    pd.(kpv{i}).tpz=po.(kpv{i}).tpz;
    pd.(kpv{i}).bpx=po.(kpv{i}).bpx;
    pd.(kpv{i}).bpy=po.(kpv{i}).bpy;
    pd.(kpv{i}).bpz=po.(kpv{i}).bpz;
    pd.(kpv{i}).cott_ft=po.(kpv{i}).lp; %cott=center of torque tube
    pd.(kpv{i}).slp_deg=po.(kpv{i}).slp;
    pd.(kpv{i}).pt_desc=string(pmd.(kpv{i}).tt);
    pmd_idx = find(contains(pd.(kpv{i}).pt_desc,"int_")==0 & ...
        contains(pd.(kpv{i}).pt_desc,"ext_")==0 & ...
        strcmp(pd.(kpv{i}).pt_desc,"NA")==0);
    pd.(kpv{i}).pt_desc(pmd_idx)=strcat("int_",string(pd.(kpv{i}).pt_desc(pmd_idx)));
    pd.(kpv{i}).pt_desc=cellstr(pd.(kpv{i}).pt_desc);
    pd.(kpv{i}).cott_cm=po.(kpv{i}).lp*30.48;
    pd.(kpv{i}).slp_pct=tand(po.(kpv{i}).slp)*100;
    pd.(kpv{i}).truss_type=po.(kpv{i}).type;
    pd.(kpv{i}).pilebin=zeros(length(pd.(kpv{i}).tpx),1);
    pd.(kpv{i}).pile_reveal=pd.(kpv{i}).cott_ft; %cott = center of torque tube height
    if const.Nevados==1
        pd.(kpv{i}).NevJoint=po.(kpv{i}).NevJoint;
        pd.(kpv{i}).LocalSlope=po.(kpv{i}).pierSumSlope;
    end
    for j=1:length(pd.(kpv{i}).pt_desc)
        if ~isempty(pd.(kpv{i}).pt_desc(j))
            if contains(pd.(kpv{i}).pt_desc(j),'mtr')==1
                pd.(kpv{i}).pile_reveal(j)=pd.(kpv{i}).pile_reveal(j)-(38/12); %32/12 from Hannah excel book 01/31/23 - should turn into input variable
            else
                pd.(kpv{i}).pile_reveal(j)=pd.(kpv{i}).pile_reveal(j)-(0.5); %top of pile is 6" below torque tube - should turn into input variable
            end
        end
    end

    pd.(kpv{i}).total_steel=zeros(length(pd.(kpv{i}).tpx),1);
    pd.(kpv{i}).xsect=cell(length(pd.(kpv{i}).tpx),1);
    [X]=discretize(pd.(kpv{i}).pile_reveal,lf.lpbin);
    for j=1:length(pd.(kpv{i}).pile_reveal)
        if ~isnan(pd.(kpv{i}).pile_reveal(j))
            pd.(kpv{i}).pilebin(j)=lf.lpbin(X(j)+1);
            pd.(kpv{i}).total_steel(j,1)=lf.lpbin(X(j)+1)+lf.embed(X(j)+1);
            pd.(kpv{i}).xsect(j,1)=lf.description(X(j)+1);
        end
    end
    pd.(kpv{i}).driven_steel=pd.(kpv{i}).total_steel-pd.(kpv{i}).pile_reveal;
    clear X

    pileTable.(kpv{i}).pile_bin=(lf.lpbin(2:length(lf.lpbin)));
    pileTable.(kpv{i}).type=(lf.type(2:length(lf.lpbin)));
    pileTable.(kpv{i}).description=(lf.description(2:length(lf.lpbin)));
    
    foo=pd.(kpv{i}).pile_reveal;
    for j=2:length(lf.lpbin)
        if contains(lf.type(j),'int')==1
            bin=foo<lf.lpbin(j);
            type=contains(pd.(kpv{i}).pt_desc,'int');
            pileTable.(kpv{i}).pile(j-1,1)=numel(foo(bin & type));
            foo(bin & type)=5000;
        elseif contains(lf.type(j),'ext')==1
            bin=foo<lf.lpbin(j);
            type=contains(pd.(kpv{i}).pt_desc,'ext');
            pileTable.(kpv{i}).pile(j-1,1)=numel(foo(bin & type));
            foo(bin & type)=5000;
        else
            pileTable.(kpv{i}).pile(j-1,1)=0;
        end
    end
    

    tln=sum(pileTable.(kpv{i}).pile(:));
    for j=1:length(pileTable.(kpv{i}).pile)
        pileTable.(kpv{i}).total_qty(j,1)=pileTable.(kpv{i}).pile(j);
        pileTable.(kpv{i}).pr(j,1)=(pileTable.(kpv{i}).total_qty(j)/tln)*100;
    end
    pileTable.(kpv{i}).pc_base=pileTable.(kpv{i}).total_qty;
    pileTable.(kpv{i})=struct2table(pileTable.(kpv{i}));
    pdout.(kpv{i})=struct2table(pd.(kpv{i}));
    pdout.(kpv{i})=rmmissing(pdout.(kpv{i}));
    oos_data.(kpv{i}).bsr=unique(pdout.(kpv{i}).bsr,'stable'); %if an error occurs here there are duplicate BSR's - check one file before loop where error occurs
    %--------------------------------------------
    oos_table=struct2table(oos_data.(kpv{i})); %if you get an error in this 
    % line make sure your truss type names are correct
    %there also appears to be an edge case where (maybe just Nevados?)
    %piles are falling outside of the surface boundary
    if i==1
        ot = oos_table;
    else
        ot = [ot;oos_table];
    end
        ot = rmmissing(ot);
    %--------------------------------------------
end
if const.writefiles==1
    for i=1:length(kpv)
        if const.simple_out==0
            mkdir([const.outpath '/output/general/' char(const.sect_name) '/' kpv{i} '_output']);
            writetable(pileTable.(kpv{i}),[const.outpath '/output/general/' char(const.sect_name) '/' kpv{i} '_output/' kpv{i} '_Pile_BOM.csv']);
        else
            mkdir([const.outpath '/' kpv{i} '_BoM/']);
            writetable(pileTable.(kpv{i}),[const.outpath '/' kpv{i} '_BoM/' kpv{i} '_Pile_BOM.csv']);
        end
    end
end
if const.writefiles==1
    for i=1:length(kpv)
        if const.simple_out==0
            mkdir([const.outpath '/output/general/' char(const.sect_name) '/' kpv{i} '_output']);
            writetable(pdout.(kpv{i}),[const.outpath '/output/general/' char(const.sect_name) '/' kpv{i} '_output/' kpv{i} '_Pile_Data_all.csv']);
        else
            mkdir([const.outpath '/' kpv{i} '_BoM/']);
            writetable(pdout.(kpv{i}),[const.outpath '/' kpv{i} '_BoM/' kpv{i} '_Pile_Data_all.csv']);
            mkdir([const.outpath '/output/.xtemp/']);
            writetable(pdout.(kpv{i}),[const.outpath '/output/.xtemp/' kpv{i} '_Pile_Data_all.csv']);
        end
    end
end
end
