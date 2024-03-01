function [drowh, dpileh]=af_2_pgplc_ojjo_6(const, surface, dpileh, drowh)
%take normal cott points
%add workpoint offset
%find leg entry based off of workpoint offset
%find pdist from workpoint to ground entry to find raw leg length
%remove workpoint to leg and reveals from raw leg length
%make work simultaneously for nxt and ati

geo=table2struct(const.geo_inputs{1});
m2f=3.28084;


%find within pdatah different types of points - contains all Ojjo point descriptions
mpts=find(contains(dpileh.pt,'m','IgnoreCase',true)); %motor points
tpts=find(~contains(dpileh.pt,'m','IgnoreCase',true)); %truss points
stpts=find(contains(dpileh.pt,'st','IgnoreCase',true)&~contains(dpileh.pt,'m','IgnoreCase',true));%standard truss
edpts=find(contains(dpileh.pt,'ed','IgnoreCase',true)&~contains(dpileh.pt,'m','IgnoreCase',true));%edge truss
hvpts=find(contains(dpileh.pt,'h','IgnoreCase',true)&~contains(dpileh.pt,'m','IgnoreCase',true));%heavy truss
hedpts=find(contains(dpileh.pt,'ed','IgnoreCase',true)&contains(dpileh.pt,'he','IgnoreCase',true));%heavy edge truss
stmpts=find(contains(dpileh.pt,'st','IgnoreCase',true)&contains(dpileh.pt,'m','IgnoreCase',true));%standard motor
hvmpts=find(contains(dpileh.pt,'h','IgnoreCase',true)&contains(dpileh.pt,'m','IgnoreCase',true));%heavy motor

%calculate the wp height above ground
lt= {mpts,                tpts}; %motor and truss points
wpv={geo.cott_wp_vert_mo, geo.cott_wp_vert_tr}; %motor and truss cott to wp
for i=1:numel(lt)
    dpileh.wpzc_ojpg(lt{i})=dpileh.tpzc(lt{i})+wpv{i}*m2f; %wp elevation z
    dpileh.wphc_ojpg(lt{i})=dpileh.cott_pg(lt{i})+wpv{i}*m2f;  %wp height
end
%calculate the y location of theoretical pile entry
for i=1:height(drowh)
    for j=drowh.si(i):drowh.ei(i)
        dpileh.bpyc_ojpg(j)=dpileh.tpyc(j)+sind(drowh.slp(i))*dpileh.wphc_ojpg(j);
    end
end

%leg insert locations xyz -- std/ed/heavy -- sttmo/hvmot
lt= {mpts,                tpts}; %motor and truss points
ang={geo.mo_angle,        geo.tr_angle}; %motor and truss angles
wpv={geo.cott_wp_vert_mo, geo.cott_wp_vert_tr}; %motor and truss cott to wp
for i=1:numel(lt)
    dpileh.xwl_ojpg(lt{i})=dpileh.tpxc(lt{i})-tand(ang{i})*dpileh.wphc_ojpg(lt{i}); %in feet
    dpileh.ywl_ojpg(lt{i})=dpileh.bpyc_ojpg(lt{i}); %legs are tangent to plane with slope of row
    dpileh.zwl_ojpg(lt{i})=surface.Fg(dpileh.xwl_ojpg(lt{i}),dpileh.ywl_ojpg(lt{i})); %surface at leg entry
    dpileh.xel_ojpg(lt{i})=dpileh.tpxc(lt{i})+tand(ang{i})*dpileh.wphc_ojpg(lt{i}); %in feet
    dpileh.yel_ojpg(lt{i})=dpileh.bpyc_ojpg(lt{i}); %legs are tangent to plane with slope of row
    dpileh.zel_ojpg(lt{i})=surface.Fg(dpileh.xel_ojpg(lt{i}),dpileh.yel_ojpg(lt{i})); %surface at leg entry
end

%find leg lengths
wp=[dpileh.tpxc,dpileh.tpyc,dpileh.wpzc_ojpg];%work point top
we=[dpileh.xwl_ojpg,dpileh.ywl_ojpg,dpileh.zwl_ojpg];%west entry point
ee=[dpileh.xel_ojpg,dpileh.yel_ojpg,dpileh.zel_ojpg];%east entry point
for i=1:height(wp)
    rwl(i,1)=pdist2(wp(i,:),we(i,:)); rel(i,1)=pdist2(wp(i,:),ee(i,:)); %raw west leg, raw east leg
end

%leg lengths - all in feet -- standard, edge, heavy, heavy edge, st motor, hv motor
%might seem like a faff to do it this way, but makes keeping things aligned easier
tt= {stpts,        edpts,        hvpts,        hedpts,        stmpts,       hvmpts}; %truss type
wpl={geo.twpl,     geo.twpl,     geo.twpl,     geo.twpl,      geo.mwpl,     geo.mwpl}; %truss wp to leg
rv= {geo.revealst, geo.revealet, geo.revealht, geo.revealhet, geo.revealsm, geo.revealhm}; %screw reveals
for i=1:numel(tt)
    dpileh.wll_ojpg(tt{i}) = rwl(tt{i})-wpl{i}*m2f-rv{i}*m2f;  %west leg length in feet
    dpileh.ell_ojpg(tt{i}) = rel(tt{i})-wpl{i}*m2f-rv{i}*m2f; %east leg length in feet
end
dpileh.wllm_ojpg=dpileh.wll_ojpg./m2f;
dpileh.ellm_ojpg=dpileh.ell_ojpg./m2f;
dpileh.pile_rev_pg=dpileh.cott_pg;
%discretize and roll up leg data
lf=const.leg_factory{1};
[X]=discretize(dpileh.wllm_ojpg,lf.lpbin);
[Y]=discretize(dpileh.ellm_ojpg,lf.lpbin);

for j=1:length(dpileh.wllm_ojpg)
    dpileh.LBW_ojpg(j)=lf.lpbin(X(j)+1);
    dpileh.LBE_ojpg(j)=lf.lpbin(Y(j)+1);
end

clear X Y
lfout.leg_bin=(lf.lpbin(2:length(lf.lpbin)));
lfout.pn=(lf.pn(2:length(lf.lpbin)));
lfout.description=(lf.description(2:length(lf.lpbin)));

for j=2:length(lf.lpbin)
    lfout.westleg(j-1,1)=numel(find(abs(dpileh.LBW_ojpg-lf.lpbin(j))<.05));
    lfout.eastleg(j-1,1)=numel(find(abs(dpileh.LBE_ojpg-lf.lpbin(j))<.05));
end

tln=sum(lfout.westleg(:))+sum(lfout.eastleg(:));

for j=1:length(lfout.westleg)
    lfout.total_qty(j,1)=lfout.westleg(j)+lfout.eastleg(j);
    lfout.lbr(j,1)=(lfout.total_qty(j)/tln)*100;
end
lfout=struct2table(lfout);
writetable(lfout,append(const.dpath{1},'/legtable_postgrade.xlsx'))
writetable(dpileh,append(const.dpath{1},'/Pile_Data_All.xlsx'))
writetable(drowh,append(const.dpath{1},'/Row_Data_All.xlsx'))

end
