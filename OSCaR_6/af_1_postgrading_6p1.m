function [drowh, dpileh] = af_1_postgrading_6p1(const, surface, drowh, dpileh)

dpileh.bpz_pg=surface.Fg(dpileh.tpxc,dpileh.tpyc); %for all foundations
if sum(strcmpi('cottcg',dpileh.Properties.VariableNames))>0
    dpileh.tpzc=dpileh.bpzc+dpileh.cottcg;
    dpileh.cott_pg=dpileh.tpzc-dpileh.bpz_pg; 
else
    dpileh.cott_pg=dpileh.tpzc-dpileh.bpz_pg;
end
others={'ATI','NXT','XTR','XTR1p5','NEV','FTC'};

if contains(char(const.tracker),'ojjo','IgnoreCase',true)
    [drowh, dpileh]=af_2_pgplc_ojjo_6(const, surface, dpileh, drowh);
elseif sum(strcmpi(char(const.tracker),others))>0
    [drowh, dpileh]=af_3_pgplc_atinxt_6(const, surface, dpileh, drowh);
end


end
