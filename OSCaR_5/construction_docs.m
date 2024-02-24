function [con_out] = construction_docs(kpv, pmd, po)
for i=1:length(kpv)
    con.(kpv{i}).bsr=pmd.(kpv{i}).bsr';
    con.(kpv{i}).id_number=pmd.(kpv{i}).id_number';
    con.(kpv{i}).section_number=pmd.(kpv{i}).sect;
    con.(kpv{i}).row_number=pmd.(kpv{i}).row;
    con.(kpv{i}).truss_letter=pmd.(kpv{i}).tl;
    con.(kpv{i}).truss_type=pmd.(kpv{i}).tt;
    con.(kpv{i}).WP_cm=po.(kpv{i}).lp*30.48;
    con.(kpv{i}).TD3=con.(kpv{i}).WP_cm>=245;
    con.(kpv{i}).slp_pct=tand(po.(kpv{i}).slp)*100;
    con.(kpv{i}).span_type=pmd.(kpv{i}).span_type;
    con_out.(kpv{i})=struct2table(con.(kpv{i}));
    con_out.(kpv{i})=rmmissing(con_out.(kpv{i}));
end