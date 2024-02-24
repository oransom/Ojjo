function [rpcs_dat,rp_loc]=rpcs_table(kpv, ot, const, pd_out, pileTable, mover_table, con_out, ppp_out, total_grad, flipex, ns_ugly, oos_master);
%clf
%pl is pile location, pls is pile location sorted
%customer/client
%project
%block/section
%row number
%pile letter
%og x
%og y
%og z - customer provided
%final x
%final y (slope adjusted)
%final z
%slope percent
%slope degree
%reveal bind
%pile reveal
%int ext
%bearing/ gear
%BHA/Motor
%std/hev/exhev
%pile size
%pile length
%size length
%upsized size
%up length
%up size length
%nxt row types
%damper y/n
%mod clear at full tilt

for i=1:length(kpv)
        for j=1:length(con_out.(kpv{i}).bsr)
            pl(i).customer_project(j,1)=const.sect_name;
        end
        pl(i).sect=pd_out.(kpv{i}).section;
        pl(i).row=pd_out.(kpv{i}).row_number;
        pl(i).pile=pd_out.(kpv{i}).truss_letter;
        pl(i).tpx_og=pd_out.(kpv{i}).tpx;
        pl(i).tpy_og=ppp_out.(kpv{i}).nsy;
        pl(i).tpz_og=pd_out.(kpv{i}).tpz;
        pl(i).bpx_og=pd_out.(kpv{i}).bpx;
        pl(i).bpy_og=ppp_out.(kpv{i}).nsy;
        pl(i).bpz_og=pd_out.(kpv{i}).bpz;
        pl(i).tpxf=ppp_out.(kpv{i}).tpx;
        pl(i).tpyf=ppp_out.(kpv{i}).tpy;
        pl(i).tpzf=ppp_out.(kpv{i}).top_pile_elev_ft;
        pl(i).bpxf=ppp_out.(kpv{i}).tpx;
        pl(i).bpyf=ppp_out.(kpv{i}).bpy;
        pl(i).bpzf=ppp_out.(kpv{i}).top_pile_elev_ft-ppp_out.(kpv{i}).reveal_ht_ft;
        pl(i).y_delta=ppp_out.(kpv{i}).bpydiff;
        pl(i).slope_pct=pd_out.(kpv{i}).slp_pct;
        pl(i).slope_deg=pd_out.(kpv{i}).slp_deg;
        pl(i).cott_ft=pd_out.(kpv{i}).cott_ft;
        %pl(i).cott_ft_frac=rat(pd_out.(kpv{i}).cott_ft,0.0625);
        pl(i).reveal_bin=pd_out.(kpv{i}).pilebin;
        pl(i).reveal=pd_out.(kpv{i}).pile_reveal;
        %pl(i).reveal_frac=rat(pd_out.(kpv{i}).pile_reveal,0.0625);
        pl(i).int_ext=pd_out.(kpv{i}).pt_desc;
        pl(i).x_sect=pd_out.(kpv{i}).xsect;
        pl(i).pile_len=pd_out.(kpv{i}).total_steel;
        for j=1:length(pl(i).x_sect)
            xsh(j)=string(pd_out.(kpv{i}).xsect(j));
            pl(i).p_size_len(j,1)=strjoin([num2str(pd_out.(kpv{i}).total_steel(j)) '_' xsh(j)]);
        end
        pl(i).up_x_sect=pd_out.(kpv{i}).upsized_xsect;
        pl(i).up_len=pd_out.(kpv{i}).total_steel;
        for j=1:length(pl(i).x_sect)
            xsh(j)=string(pd_out.(kpv{i}).upsized_xsect(j));
            pl(i).up_p_size_len(j,1)=strjoin([num2str(pd_out.(kpv{i}).total_steel(j)) '_' xsh(j)]);
        end
        if i==1
            pl=pl(i);
        else
            pl=[pl,pl(i)];
        end
end

rpcs_dat=pl;
rpcs_dat_out=struct2table(pl);
if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/' kpv{i} '_output']);
        writetable(rpcs_dat_out,[const.outpath '/' kpv{i} '_output/' kpv{i} '_RPCS_DATA.csv'])
        rp_loc=[const.outpath '/' kpv{i} '_output/' kpv{i} '_RPCS_DATA.csv'];
    else
        mkdir([const.outpath '/' char(const.project{1}) '_output_' char(const.t{1}) '/rpdf']);
        writetable(rpcs_dat_out,[const.outpath '/' char(const.project{1}) '_output_' char(const.t{1}) '/rpdf/' char(const.sect_name) '_RPCS_DATA.csv'])
        rp_loc=[const.outpath '/' char(const.project{1}) '_output_' char(const.t{1}) '/rpdf/' char(const.sect_name) '_RPCS_DATA.csv'];
    end
end