function [pd_out,piletable_upr,ppp_out,S_sI_pg] = postdata_postgrade(kpv, graded_xyz, pd_out, const, ppp_out)
S_sI_pg=scatteredInterpolant(graded_xyz(:,1),graded_xyz(:,2),graded_xyz(:,3),'natural','none');
lf=const.leg_factory{1};
pb.input=const.pile_bins{1};
pb.input=table2cell(pb.input);
mtr_wp2top=const.mtr_wp2top; %32/12 from Hannah excel book 01/31/23 - should turn into input variable
std_wp2top=const.std_wp2top; %top of pile is 6" below torque tube - should turn into input variable
%spg=S_sI_pg(surface.xq,surface.yq);
for i=1:length(kpv)
    pd.(kpv{i})=pd_out.(kpv{i});
    pd.(kpv{i}).bpy=pd.(kpv{i}).tpy;
    pd.(kpv{i}).bpz=S_sI_pg(pd.(kpv{i}).bpx,pd.(kpv{i}).bpy);
    pd.(kpv{i}).cott_ft=pd.(kpv{i}).tpz-pd.(kpv{i}).bpz;
    pd.(kpv{i}).cott_cm=pd.(kpv{i}).cott_ft*30.48;
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
        upsized_xsect_array=upsized_xsect_array(~cellfun('isempty',upsized_xsect_array));
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
                %comment out next line for ojjo
                %pd.(kpv{i}).upsized_pilebin(j)=lf.lpbin(min(newbin_ind(ind))); %you're getting an error here because the span names aren't 7 chars long with the first 3 being int or ex
                % or your different spans aren't restarting numbering at 1
                pd.(kpv{i}).upsized_pilebin(j)=pd.(kpv{i}).upsized_pilebin(j);
            end
            %end
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

    pd_out.(kpv{i})=pd.(kpv{i});
    pd_out.(kpv{i})=rmmissing(pd_out.(kpv{i}));

    % %% NEW PPP
    ppp_out.(kpv{i}).tpx=pd.(kpv{i}).tpx;
    ppp_out.(kpv{i}).tpy=pd.(kpv{i}).tpy;
    ppp_out.(kpv{i}).reveal_ht_ft=pd.(kpv{i}).cott_ft;
    ppp_out.(kpv{i}).reveal_ht_cm=pd.(kpv{i}).cott_cm;
    ppp_out.(kpv{i}).top_pile_elev_ft=pd.(kpv{i}).tpz;
    ppp_out.(kpv{i}).bpy=pd.(kpv{i}).bpy;
    ppp_out.(kpv{i}).bpydiff=pd.(kpv{i}).bpy-ppp_out.(kpv{i}).nsy;
    ppp_out.(kpv{i}).slope_pct=pd.(kpv{i}).slp_pct;

    if const.writefiles==1
        mkdir([const.outpath '/output/']);
        for i=1:length(kpv)
            writetable(pd_out.(kpv{i}),[const.outpath '/output/' 'Pile_Data_all.csv']);
        end
    end
end
end
