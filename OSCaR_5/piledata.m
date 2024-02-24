function [pdout, lfout, ot, piletable_upr] = piledata(kpv, pmd, po, to, oos_data, const)
if contains(const.tracker,'ojjo','IgnoreCase',true)
    geo=table2struct(const.geo_inputs{1}); %going to have to do something with this for non ojjo 12/11/23 - maybe I did w/ if statement?
end
lf=const.leg_factory{1};
pb.input=const.pile_bins{1};
pb.input=table2cell(pb.input);
mtr_wp2top=const.mtr_wp2top; %32/12 from Hannah excel book 01/31/23 - should turn into input variable
std_wp2top=const.std_wp2top; %top of pile is 6" below torque tube - should turn into input variable
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
    if contains(const.tracker,'ojjo','IgnoreCase',true)
        pd.(kpv{i}).L_W_ft=to.(kpv{i}).lw;
        pd.(kpv{i}).WP_ft=po.(kpv{i}).lp;
        pd.(kpv{i}).L_E_ft=to.(kpv{i}).le;
        if const.thirdleg==1
            pd.(kpv{i}).TL_ft=to.(kpv{i}).lm;
        end
        pd.(kpv{i}).L_W_m=to.(kpv{i}).lw*.3048;
        pd.(kpv{i}).WP_cm=po.(kpv{i}).lp*30.48;
        pd.(kpv{i}).L_E_m=to.(kpv{i}).le*.3048;
        if const.thirdleg==1
            pd.(kpv{i}).TL_m=to.(kpv{i}).lm*.3048;
        end
    end
    pd.(kpv{i}).cott_ft=po.(kpv{i}).lp; %cott=center of torque tube
    pd.(kpv{i}).slp_deg=po.(kpv{i}).slp;
    pd.(kpv{i}).pt_desc=string(pmd.(kpv{i}).tt);
    if ~contains(const.tracker,'ojjo','IgnoreCase',true)
        pmd_idx = find(contains(pd.(kpv{i}).pt_desc,"int_")==0 & ...
            contains(pd.(kpv{i}).pt_desc,"ext_")==0 & ...
            strcmp(pd.(kpv{i}).pt_desc,"NA")==0);
        pd.(kpv{i}).pt_desc(pmd_idx)=strcat("int_",string(pd.(kpv{i}).pt_desc(pmd_idx)));
        pd.(kpv{i}).pt_desc=cellstr(pd.(kpv{i}).pt_desc);
    end
    pd.(kpv{i}).cott_cm=po.(kpv{i}).lp*30.48;
    pd.(kpv{i}).slp_pct=tand(po.(kpv{i}).slp)*100;
    pd.(kpv{i}).truss_type=po.(kpv{i}).type;
    if contains(const.tracker,'ojjo','IgnoreCase',true)
        pd.(kpv{i}).ALW=zeros(length(pd.(kpv{i}).truss_type),1);
        pd.(kpv{i}).LBW=zeros(length(pd.(kpv{i}).truss_type),1);
        pd.(kpv{i}).ALE=zeros(length(pd.(kpv{i}).truss_type),1);
        pd.(kpv{i}).LBE=zeros(length(pd.(kpv{i}).truss_type),1);
        if const.thirdleg==1
            pd.(kpv{i}).ALM=zeros(length(pd.(kpv{i}).truss_type),1);
            pd.(kpv{i}).LBM=zeros(length(pd.(kpv{i}).truss_type),1);
        end
        pdh.es=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Standard"))-(geo.twpl+geo.revealst); %add it all up, in the coordinate system of the leg at an angle
        pdh.esm=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Standard Motor"))-(geo.mwpl+geo.revealsm);
        pdh.eh=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy"))-(geo.twpl+geo.revealht);
        pdh.ehm=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy Motor"))-(geo.mwpl+geo.revealhm);
        pdh.ee=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Edge"))-(geo.twpl+geo.revealet);
        pdh.eem=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Edge Motor"))-(geo.mwpl+geo.revealhm);
        pdh.ehe=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy Edge"))-(geo.mwpl+geo.revealhet);
        pdh.ws=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Standard"))-(geo.twpl+geo.revealst);
        pdh.wsm=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Standard Motor"))-(geo.mwpl+geo.revealsm);
        pdh.wh=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy"))-(geo.twpl+geo.revealht);
        pdh.whm=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy Motor"))-(geo.mwpl+geo.revealhm);
        pdh.we=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Edge"))-(geo.twpl+geo.revealet);
        pdh.wem=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Edge Motor"))-(geo.mwpl+geo.revealhm);
        pdh.whe=pd.(kpv{i}).L_W_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy Edge"))-(geo.mwpl+geo.revealhet);

        [lf_hev_mtr_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Heavy Motor")); %find indicies of locations
        [lf_hev_std_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Heavy"));
        [lf_std_mtr_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Standard Motor"));
        [lf_std_std_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Standard"));
        [lf_edg_mtr_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Edge Motor"));
        [lf_edg_std_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Edge"));
        [lf_edg_hev_idx, ~] = find(matches(string(pd.(kpv{i}).pt_desc),"Heavy Edge"));

        for j=1:length(lf_hev_mtr_idx)
            pd.(kpv{i}).ALE(lf_hev_mtr_idx(j),1)=pdh.ehm(j); %insert locations in their proper locations
            pd.(kpv{i}).ALW(lf_hev_mtr_idx(j),1)=pdh.whm(j);
        end
        for j=1:length(lf_hev_std_idx)
            pd.(kpv{i}).ALE(lf_hev_std_idx(j),1)=pdh.eh(j);
            pd.(kpv{i}).ALW(lf_hev_std_idx(j),1)=pdh.wh(j);
        end
        for j=1:length(lf_std_mtr_idx)
            pd.(kpv{i}).ALE(lf_std_mtr_idx(j),1)=pdh.esm(j);
            pd.(kpv{i}).ALW(lf_std_mtr_idx(j),1)=pdh.wsm(j);
        end
        for j=1:length(lf_std_std_idx)
            pd.(kpv{i}).ALE(lf_std_std_idx(j),1)=pdh.es(j);
            pd.(kpv{i}).ALW(lf_std_std_idx(j),1)=pdh.ws(j);
        end
        for j=1:length(lf_edg_mtr_idx)
            pd.(kpv{i}).ALE(lf_edg_mtr_idx(j),1)=pdh.eem(j);
            pd.(kpv{i}).ALW(lf_edg_mtr_idx(j),1)=pdh.wem(j);
        end
        for j=1:length(lf_edg_std_idx)
            pd.(kpv{i}).ALE(lf_edg_std_idx(j),1)=pdh.ee(j);
            pd.(kpv{i}).ALW(lf_edg_std_idx(j),1)=pdh.we(j);
        end
        for j=1:length(lf_edg_hev_idx)
            pd.(kpv{i}).ALE(lf_edg_hev_idx(j),1)=pdh.ehe(j);
            pd.(kpv{i}).ALW(lf_edg_hev_idx(j),1)=pdh.whe(j);
        end

        [X]=discretize(pd.(kpv{i}).ALW,lf.lpbin);
        [Y]=discretize(pd.(kpv{i}).ALE,lf.lpbin);
        if const.thirdleg==1
            [Z]=discretize(pd.(kpv{i}).ALM,lf.lpbin);
        end

        for j=1:length(pd.(kpv{i}).ALW)
            if ~isnan(pd.(kpv{i}).ALW(j))
                pd.(kpv{i}).LBW(j)=lf.lpbin(X(j)+1);
                pd.(kpv{i}).LBE(j)=lf.lpbin(Y(j)+1);
                if const.thirdleg==1
                    if pd.(kpv{i}).ALM(j)~=0
                        pd.(kpv{i}).LBM(j)=lf.lpbin(Z(j)+1);
                    else
                        pd.(kpv{i}).LBM(j)=0;
                    end
                end
            end
        end
        clear X Y Z
        lfout.(kpv{i}).leg_bin=(lf.lpbin(2:length(lf.lpbin)));
        lfout.(kpv{i}).pn=(lf.pn(2:length(lf.lpbin)));
        lfout.(kpv{i}).description=(lf.description(2:length(lf.lpbin)));

        for j=2:length(lf.lpbin)
            lfout.(kpv{i}).westleg(j-1,1)=numel(find(abs(pd.(kpv{i}).LBW-lf.lpbin(j))<.05));
            lfout.(kpv{i}).eastleg(j-1,1)=numel(find(abs(pd.(kpv{i}).LBE-lf.lpbin(j))<.05));
            if const.thirdleg==1
                lfout.(kpv{i}).motorleg(j-1,1)=numel(find(abs(pd.(kpv{i}).LBM-lf.lpbin(j))<0.05));
            end
        end
        if const.thirdleg==1
            tln=sum(lfout.(kpv{i}).westleg(:))+sum(lfout.(kpv{i}).eastleg(:))+sum(lfout.(kpv{i}).motorleg(:));
        else
            tln=sum(lfout.(kpv{i}).westleg(:))+sum(lfout.(kpv{i}).eastleg(:));
        end
        for j=1:length(lfout.(kpv{i}).westleg)
            if const.thirdleg==1
                lfout.(kpv{i}).total_qty(j,1)=lfout.(kpv{i}).westleg(j)+lfout.(kpv{i}).eastleg(j)+lfout.(kpv{i}).motorleg(j);
            else
                lfout.(kpv{i}).total_qty(j,1)=lfout.(kpv{i}).westleg(j)+lfout.(kpv{i}).eastleg(j);
            end
            lfout.(kpv{i}).lbr(j,1)=(lfout.(kpv{i}).total_qty(j)/tln)*100;
        end
        lfout.(kpv{i})=struct2table(lfout.(kpv{i}));
    else
        lfout.(kpv{i})=[];
    end
    pd.(kpv{i}).pilebin=zeros(length(pd.(kpv{i}).tpx),1);
    pd.(kpv{i}).pile_reveal=pd.(kpv{i}).cott_ft; %cott = center of torque tube height
    for j=1:length(pd.(kpv{i}).pt_desc)
        if ~isempty(pd.(kpv{i}).pt_desc(j))
            if contains(pd.(kpv{i}).pt_desc(j),'mtr')==1
                pd.(kpv{i}).pile_reveal(j)=pd.(kpv{i}).pile_reveal(j)-mtr_wp2top;
            else
                pd.(kpv{i}).pile_reveal(j)=pd.(kpv{i}).pile_reveal(j)-std_wp2top;
            end
        end
    end

    pd.(kpv{i}).total_steel=zeros(length(pd.(kpv{i}).tpx),1);
    pd.(kpv{i}).xsect=cell(length(pd.(kpv{i}).tpx),1);
    [X]=discretize(pd.(kpv{i}).pile_reveal,lf.lpbin); %changed this from cott to pile reveal on 12/19/23
    for j=1:length(pd.(kpv{i}).pile_reveal)
        if ~isnan(pd.(kpv{i}).pile_reveal(j))
            pd.(kpv{i}).pilebin(j)=lf.lpbin(X(j)+1);
            pd.(kpv{i}).total_steel(j,1)=lf.lpbin(X(j)+1)+lf.embed(X(j)+1);
        end
    end
    pd.(kpv{i}).driven_steel=pd.(kpv{i}).total_steel-pd.(kpv{i}).pile_reveal;
    clear X
    for j=1:length(lf.lpbin)
        mask=find(pd.(kpv{i}).pilebin==lf.lpbin(j) & contains(pd.(kpv{i}).pt_desc,lf.type{j}));
        pd.(kpv{i}).xsect(mask,1)=lf.description(j); %I can see an error happening here in the future
    end



    %rolling up pile
    if ~contains(const.tracker,'ojjo','IgnoreCase',true)
        numsect=1:max(pd.(kpv{i}).section);
        l=0;
        for j=1:max(numsect)
            numrow=unique(pd.(kpv{i}).row_number(pd.(kpv{i}).section==j));
            for k=1:length(numrow)
                kk=numrow(k);
                l=l+1;
                pilet{:,l}=pd.(kpv{i}).xsect(pd.(kpv{i}).section==j & pd.(kpv{i}).row_number==kk); %adding xsect into pile data
                pileb{:,l}=pd.(kpv{i}).pilebin(pd.(kpv{i}).section==j & pd.(kpv{i}).row_number==kk); %adding pile bin into pile data
                foosect{:,l}=pd.(kpv{i}).section(pd.(kpv{i}).section==j & pd.(kpv{i}).row_number==kk);
                foorow{:,l}=pd.(kpv{i}).row_number(pd.(kpv{i}).section==j & pd.(kpv{i}).row_number==kk);
            end
        end

        for j=1:length(pilet)
            pileuph=pilet{j}; %turn the cells to chars
            for k=1:size(pileuph,1)
                pileup=strtrim(pileuph(k,:));%took the char command out of here - hope it doesn't break anything
                pileup=char(pileup); %if you get errors here make sure your span names include (int/ext/mtr, etc.) or pile falling off surface
                [pb.ri(j,k) pb.ci(j,k)] = find(strcmp(pileup,pb.input));
            end
        end
        pb.ri(pb.ri==0)=NaN;
        pb.ci(pb.ci==0)=NaN;
        [j k] = size(pb.ri);
        for m=1:j
            maxr=max(pb.ri(m,:)); %max row value in row of pile
            minc=(pb.ci(m,:)>0 & pb.ri(m,:)==maxr).*pb.ci(m,:); %min column in that max row
            minc(minc==0)=NaN; %change zeros to nan
            mincf=min(minc); %min non-zero column in the max row
            maxcf=max(minc); %max non-zero column in the max row
            colsame=find(pb.ri(m,:)==max(pb.ri(m,:))); %columns where row = maxrow
            coldiff=find(pb.ri(m,:)<max(pb.ri(m,:))); %columns where row is less than maxrow
            rowsame=find(pb.ri(m,:)==max(pb.ri(m,:))); %rows matching max row
            rowdiff=find(pb.ri(m,:)<max(pb.ri(m,:))); %rows less than max row
            pb.rf(m,rowsame)=pb.ri(m,rowsame); %don't do anything to max rows
            pb.rf(m,rowdiff)=maxr; %make rows less than max row, max row
            pb.cf(m,colsame)=pb.ci(m,colsame); %don't change columns where row = maxrow
            pb.cf(m,coldiff)=mincf; %make column smallest column in max row if row is less than max row
            pb.cf(pb.cf==0)=NaN;
            pb.rf(pb.rf==0)=NaN;
        end
        for m=1:j
            for n=1:k
                if ~isnan(pb.rf(m,n)) && ~isnan(pb.cf(m,n))
                    upsized_xsect(m,n)=pb.input(pb.rf(m,n),pb.cf(m,n));
                end
            end
        end

        upsized_xsect_array=reshape(upsized_xsect',[],1);
        ind=find(~cellfun(@isempty,pd.(kpv{i}).xsect));
        pd.(kpv{i}).upsized_xsect=cell(length(pd.(kpv{i}).xsect),1);
        pd.(kpv{i}).upsized_xsect(ind)=upsized_xsect_array(ind);
        pd.(kpv{i}).upsized_pilebin=zeros(size(pd.(kpv{i}).pilebin,1),1);
        for j=1:length(pd.(kpv{i}).upsized_xsect)
            %if ~strcmp(pd.(kpv{i}).pt_desc(j),"NA")
            pd_hold=char(pd.(kpv{i}).pt_desc(j));
            pd_hold=pd_hold(1:end-4); %cheap way to isolate int & ext
            oldbin_ind=find(contains(lf.type,pd_hold) & (lf.lpbin==pd.(kpv{i}).pilebin(j))); %original pilebin (reveal length)=original
            newbin_ind=find(contains(lf.type,pd_hold) & strcmp(lf.description,pd.(kpv{i}).upsized_xsect(j))); %find new xsect reveal lengths (pilebin)
            if ~isempty(newbin_ind)
                ind=find(newbin_ind-oldbin_ind>=0); %find overlapping new and old bins, or minimum new bin above original bin
                %if ~contains(const.tracker,'ojjo','IgnoreCase',true)
                    %pd.(kpv{i}).upsized_pilebin(j)=lf.lpbin(min(newbin_ind(ind))); %you're getting an error here because the span names aren't 7 chars long with the first 3 being int or ex
                    % or your different spans aren't restarting numbering at 1
                %else
                    pd.(kpv{i}).upsized_pilebin(j)=pd.(kpv{i}).upsized_pilebin(j);
               % end
            end
        end

        pileTable.(kpv{i}).pile_bin=(lf.lpbin(2:length(lf.lpbin)));
        pileTable.(kpv{i}).type=(lf.type(2:length(lf.lpbin)));
        pileTable.(kpv{i}).description=(lf.description(2:length(lf.lpbin)));

        %% pile table upsized and consolidated
        foo=pd.(kpv{i}).pile_reveal+0.1;
        for j=1:size(lf.lpbin,1)
            ncount(j,1)=sum(lf.lpbin(j)==pd.(kpv{i}).pilebin & contains(pd.(kpv{i}).pt_desc,lf.type{j}) & strcmp(pd.(kpv{i}).xsect,lf.description{j}));
            if ncount(j)>0 && mod(ncount(j)/15,1)~=0
                ncount_bin(j,1)=(ncount(j)+(15-15*(mod(ncount(j)/15,1)))); %rounding to the nearest 15
            elseif mod(ncount(j)/15,1)==0
                ncount_bin(j,1)=ncount(j);
            else
                ncount_bin(j,1)=0;
            end
            upcount(j,1)=sum(lf.lpbin(j)==pd.(kpv{i}).upsized_pilebin & contains(pd.(kpv{i}).pt_desc,lf.type{j}) & strcmp(pd.(kpv{i}).upsized_xsect,lf.description{j}));
            if upcount(j,1)>0 && mod(upcount(j)/15,1)~=0
                upcount_bin(j)=(upcount(j)+(15-15*(mod(upcount(j)/15,1)))); %rounding to the nearest 15
            elseif mod(upcount(j)/15,1)==0
                upcount_bin(j)=upcount(j);
            else
                upcount_bin(j,1)=0;
            end
        end

        piletable_upr.(kpv{i}).pile_bin=lf.lpbin(2:end);
        piletable_upr.(kpv{i}).type=lf.type(2:end);
        piletable_upr.(kpv{i}).description=lf.description(2:end);
        piletable_upr.(kpv{i}).pc_base=ncount(2:end);
        piletable_upr.(kpv{i}).pc_base_mod=ncount_bin(2:end);
        piletable_upr.(kpv{i}).pc_consol=upcount(2:end);
        piletable_upr.(kpv{i}).pc_consol_mod=(upcount_bin(2:end))';
        piletable_upr.(kpv{i}).pc_base_pct=(ncount(2:end)/sum(ncount(2:end)));
        piletable_upr.(kpv{i}).pc_consol_pct=(upcount(2:end)/sum(upcount(2:end)));
        piletable_upr.(kpv{i}).pc_consol_mod_pct=(upcount_bin(2:end)/sum(upcount_bin(2:end)))';
        piletable_upr.(kpv{i})=struct2table(piletable_upr.(kpv{i}));
    else
        piletable_upr=[];
    end
    pdout.(kpv{i})=struct2table(pd.(kpv{i}));
    pdout.(kpv{i})=rmmissing(pdout.(kpv{i}));
    oos_data.(kpv{i}).bsr=unique(pdout.(kpv{i}).bsr,'stable'); %if an error occurs here there are duplicate BSR's - check one file before loop where error occurs
    %--------------------------------------------
    oos_table=struct2table(oos_data.(kpv{i})); %if you get an error in this line make sure your truss type names are correct
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
        if ~contains(const.tracker,'ojjo','IgnoreCase',true)
            mkdir([const.outpath '/ojjo/' kpv{i} '_BoM']);
            writetable(piletable_upr.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_BoM/' kpv{i} '_Pile_BOM.csv']);
        else
            mkdir([const.outpath '/output/'])
            writetable(pdout.(kpv{i}),[const.outpath '/output/' kpv{i} '_Pile_Data_all.csv']);
        end
    end
end
if contains(const.tracker,'ojjo','IgnoreCase',true)
    if const.writefiles==1
        for i=1:length(kpv)
            mkdir([const.outpath '/ojjo/' kpv{i} '_LegFact']);
            writetable(lfout.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_LegFact/' kpv{i} '_legFactory_pregrade.csv']);
        end
    end
end
end