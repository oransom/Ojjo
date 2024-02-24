function [pd_out,con_out,ppp_out,S_sI_pg,top_out,pin_out] = legdata_postgrade(kpv, graded_xyz, pd_out, const, con_out, ppp_out, top_out, pin_out)
if contains(const.tracker,'ojjo','IgnoreCase',true)
    geo=table2struct(const.geo_inputs{1}); %going to have to do something with this for non ojjo 12/11/23 - maybe I did w/ if statement?
end
lf=const.leg_factory{1};
cm_ft=2.54*12; %for converting
S_sI_pg=scatteredInterpolant(graded_xyz(:,1),graded_xyz(:,2),graded_xyz(:,3),'natural','none');
lf=const.leg_factory{1};
pb.input=const.pile_bins{1};
pb.input=table2cell(pb.input);
mtr_wp2top=const.mtr_wp2top; %32/12 from Hannah excel book 01/31/23 - should turn into input variable
std_wp2top=const.std_wp2top; %top of pile is 6" below torque tube - should turn into input variable
%spg=S_sI_pg(surface.xq,surface.yq);
for i=1:length(kpv)
    pd.(kpv{i})=pd_out.(kpv{i});
    pdh.zadj=ppp_out.(kpv{i}).reveal_ht_ft-pd.(kpv{i}).WP_ft;
    pd.(kpv{i}).WP_ft=pd.(kpv{i}).tpz-(S_sI_pg(pd.(kpv{i}).tpx,pd.(kpv{i}).tpy));%+geo.cott_wp_vert*3.281;
    for j=1:length(pd.(kpv{i}).tpz)
        if contains(pd.(kpv{i}).pt_desc(j),'motor','IgnoreCase',true)
            pd.(kpv{i}).tpz(j)=pd.(kpv{i}).tpz(j)+pdh.zadj(j)+(geo.cott_wp_vert_mo*3.281);
        else
            pd.(kpv{i}).tpz(j)=pd.(kpv{i}).tpz(j)+pdh.zadj(j)+(geo.cott_wp_vert_tr*3.281);
        end
    end
    pd.(kpv{i}).WP_ft=pd.(kpv{i}).tpz-(S_sI_pg(pd.(kpv{i}).tpx,pd.(kpv{i}).tpy));%+geo.cott_wp_vert*3.281;

    if const.trusscalc==1
        pd.(kpv{i}).bpy=pd.(kpv{i}).tpy+sind(pd.(kpv{i}).slp_deg).*(pd.(kpv{i}).WP_ft);
    else
        pd.(kpv{i}).bpy=pd.(kpv{i}).tpy;
    end
    pd.(kpv{i}).bpz=S_sI_pg(pd.(kpv{i}).bpx,pd.(kpv{i}).bpy);
    if contains(const.tracker,'ojjo','IgnoreCase',true)
        pdh.ex=pd.(kpv{i}).tpx+tand(const.truss_leg_angle)*pd.(kpv{i}).WP_ft;
        pdh.ey=pd.(kpv{i}).bpy;
        pdh.ez=S_sI_pg(pdh.ex,pdh.ey);
        pdh.wx=pd.(kpv{i}).tpx-tand(const.truss_leg_angle)*pd.(kpv{i}).WP_ft;
        pdh.wy=pd.(kpv{i}).bpy;
        pdh.wz=S_sI_pg(pdh.wx,pdh.wy);
    end

    pdh.t=[pd.(kpv{i}).tpx,pd.(kpv{i}).tpy,pd.(kpv{i}).tpz];
    if contains(const.tracker,'ojjo','IgnoreCase',true)
        pdh.e=[pdh.ex,pdh.ey,pdh.ez];
        pdh.w=[pdh.wx,pdh.wy,pdh.wz];

        for j=1:size(pdh.t,1)
            pd.(kpv{i}).L_W_ft(j)=pdist2(pdh.t(j,:),pdh.w(j,:));
            pd.(kpv{i}).L_W_m(j)=pdist2(pdh.t(j,:),pdh.w(j,:))*cm_ft/100;
            pd.(kpv{i}).L_E_ft(j)=pdist2(pdh.t(j,:),pdh.e(j,:));
            pd.(kpv{i}).L_E_m(j)=pdist2(pdh.t(j,:),pdh.e(j,:))*cm_ft/100;
        end
    end
    if contains(const.tracker,'ojjo','IgnoreCase',true)
        pd.(kpv{i}).WP_cm=pd.(kpv{i}).WP_ft*cm_ft;
        pdh.es=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Standard"))-(geo.twpl+geo.revealst);
        pdh.esm=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Standard Motor"))-(geo.mwpl+geo.revealsm);
        pdh.eh=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy"))-(geo.twpl+geo.revealht);
        pdh.ehm=pd.(kpv{i}).L_E_m(matches(string(pd.(kpv{i}).pt_desc),"Heavy Motor"))-(geo.mwpl+geo.revealhm); %add it all up, in the coordinate system of the leg at an angle
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

        clear X Y Z ldh

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
    pd_out.(kpv{i})=pd.(kpv{i});
    pd_out.(kpv{i})=rmmissing(pd_out.(kpv{i}));

    %% NEW PPP
    ppp_out.(kpv{i}).tpx=pd.(kpv{i}).tpx;
    ppp_out.(kpv{i}).tpy=pd.(kpv{i}).tpy;
    ppp_out.(kpv{i}).reveal_ht_ft=pd.(kpv{i}).WP_ft;
    ppp_out.(kpv{i}).reveal_ht_cm=pd.(kpv{i}).WP_cm;
    ppp_out.(kpv{i}).top_pile_elev_ft=pd.(kpv{i}).tpz;
    ppp_out.(kpv{i}).bpy=pd.(kpv{i}).bpy;
    ppp_out.(kpv{i}).bpydiff=pd.(kpv{i}).bpy-ppp_out.(kpv{i}).nsy;
    ppp_out.(kpv{i}).slope_pct=pd.(kpv{i}).slp_pct;
    %%% NEW CONRAW
    con_out.(kpv{i}).WP_cm=pd.(kpv{i}).WP_cm;
    con_out.(kpv{i}).WP_cm_final=pd.(kpv{i}).WP_cm;
    %% NEW TOPOUT
    top_out.(kpv{i}).n=pd.(kpv{i}).tpy;
    top_out.(kpv{i}).e=pd.(kpv{i}).tpx;
    top_out.(kpv{i}).z=pd.(kpv{i}).tpz;
    %% NEW PINOUT
    pin_out.(kpv{i}).n=pd.(kpv{i}).bpy;
    pin_out.(kpv{i}).e=pd.(kpv{i}).bpx;
    pin_out.(kpv{i}).z=pd.(kpv{i}).bpz;
end
if const.trusscalc==1
    if const.writefiles==1
        for i=1:length(kpv)
            mkdir([const.outpath '/ojjo/' kpv{i} '_output/']);
            writetable(lfout.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_output/' '1_LegFactory_postgrade.csv']);
            writetable(pd_out.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_output/' '1_PG_Pile_Data_all.csv']);
            writetable(ppp_out.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_output/' '2_pier-plot-plan-raw.csv']);
            writetable(con_out.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_output/' '2_con-raw.csv']);
            writetable(top_out.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_output/' '2_top_pile_xyz.csv']);
            writetable(pin_out.(kpv{i}),[const.outpath '/ojjo/' kpv{i} '_output/' '2_pin-master.csv']);
        end
    end
end
end
