function [master_dat]=mastercombiner_pile(kpv, ot, const, pd_out, lfout, pileTable, mover_table, con_out, ppp_out, total_grad, flipex, ns_ugly, oos_master, grd_area);
%clf
%pl is pile location, pls is pile location sorted
for i=1:length(kpv)
    if i==1
        pl.bsr=con_out.(kpv{i}).bsr;
        pl.tpx=pd_out.(kpv{i}).tpx;
        pl.tpy=pd_out.(kpv{i}).tpy;
        pl.tpz=pd_out.(kpv{i}).tpz;
        pl.bpx=pd_out.(kpv{i}).bpx;
        pl.bpy=pd_out.(kpv{i}).bpy;
        pl.bpz=pd_out.(kpv{i}).bpz;
        if ~contains(const.tracker,'ojjo','IgnoreCase',true)
        pl.pilecount=sum(pileTable.(kpv{i}).pc_base);
        end
        pl.steelQ=sum(pd_out.(kpv{i}).total_steel);
        if contains(const.tracker,'ojjo','IgnoreCase',true)
            if const.trusscalc==0
                pl.pilecount=sum(lfout.(kpv{i}).westleg);
            elseif const.trusscalc==1
                if const.thirdleg==1
                    pl.motorlegcount=sum(lfout.(kpv{i}).motorleg);
                end
                pl.trusscount=sum(lfout.(kpv{i}).total_qty);
                pl.trusssteel=sum(pd_out.(kpv{i}).L_W_ft)+sum(pd_out.(kpv{i}).L_E_ft);
            end
        end
    else
        pl.bsr=[pl.bsr;con_out.(kpv{i}).bsr];
        pl.tpx=[pl.tpx;pd_out.(kpv{i}).tpx];
        pl.tpy=[pl.tpy;pd_out.(kpv{i}).tpy];
        pl.tpz=[pl.tpz;pd_out.(kpv{i}).tpz];
        pl.bpx=[pl.bpx;pd_out.(kpv{i}).bpx];
        pl.bpy=[pl.bpy;pd_out.(kpv{i}).bpy];
        pl.bpz=[pl.bpz;pd_out.(kpv{i}).bpz];
        if ~contains(const.tracker,'ojjo','IgnoreCase',true)
        pl.pilecount=pl.pilecount+sum(pileTable.(kpv{i}).pc_base);
        end
        pl.steelQ=pl.steelQ+sum(pd_out.(kpv{i}).total_steel);
        if contains(const.tracker,'ojjo','IgnoreCase',true)
            if const.trusscalc==0
                pl.pilecount=pl.pilecount+sum(lfout.(kpv{i}).westleg);
            elseif const.trusscalc==1
                if const.thirdleg==1
                    pl.motorlegcount=pl.motorlegcount+sum(lfout.(kpv{i}).motorleg);
                end
                pl.trusscount=pl.trusscount+sum(lfout.(kpv{i}).total_qty);
                pl.trusssteel=pl.trusssteel+sum(pd_out.(kpv{i}).L_W_ft)+sum(pd_out.(kpv{i}).L_E_ft);
            end
        end
    end
end

for i=1:length(kpv)
      for j=1:length(pl.bsr)
          if find(strcmp(pl.bsr(j),mover_table.bsr))>0
              foo = find(strcmp(pl.bsr(j),mover_table.bsr));
              pl.tpzn(j,1)=pl.tpz(j)+mover_table.raiserow(foo);
          else
              pl.tpzn(j,1)=pl.tpz(j);
          end
      end
end

pl.tpz=pl.tpzn;    clear pl_fin.tpzn foo;

master_dat.run=char(const.sect_name);
if contains(const.tracker,'ojjo','IgnoreCase',true)
    if const.trusscalc==0
        master_dat.total_ag_pile_steel=sum(pl.tpz-pl.bpz);
        master_dat.pilecount=pl.pilecount;
    elseif const.trusscalc==1
        master_dat.trusscount=pl.trusscount;
        if const.thirdleg==1
            master_dat.motorlegcount=pl.motorlegcount;
        end
        master_dat.total_ag_truss_steel=pl.trusssteel;
    end
end
master_dat.total_ag_pile_steel=sum(pl.tpz-pl.bpz);
if ~contains(const.tracker,'ojjo','IgnoreCase',true)
master_dat.pilecount=pl.pilecount;
end
master_dat.total_steel=pl.steelQ;
master_dat.cut_cy=total_grad.cut;
master_dat.fill_cy=total_grad.fill;
master_dat.grd_area=grd_area;

if isempty(flipex.flip)==0
    master_dat.flipex=sum(flipex.flip);
else
    master_dat.flipex=0;
end
if sum(oos_master(:,6))>1
    master_dat.pilerequiringgrading=nnz(oos_master(:,6));
end
master_dat.rowssolvedwithsteel=sum(ot.preveal_max<const.max_wp & ot.preveal_min>const.min_wp);
if sum(ns_ugly.nsd)~=0
    master_dat.northsouthugly=nnz(ns_ugly.nsd);
end
master_dat_out=struct2table(master_dat);
if const.writefiles == 1
    mkdir([const.outpath '/rollups']);
    writetable(master_dat_out,[const.outpath '/rollups/' char(const.sect_name) '_DRU.csv'])
    mkdir([const.outpath '/output/.xtemp/']);
    writetable(master_dat_out,[const.outpath '/output/.xtemp/' char(const.sect_name) '_DRU.csv'])
end