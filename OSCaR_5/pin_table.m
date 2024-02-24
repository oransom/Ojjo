function [pin_out] = pin_table(kpv, pmd, po)

for i=1:length(kpv)
    pin.(kpv{i}).p=pmd.(kpv{i}).bsrt';
    pin.(kpv{i}).n=po.(kpv{i}).bpy;
    pin.(kpv{i}).e=po.(kpv{i}).bpx;
    pin.(kpv{i}).z=po.(kpv{i}).bpz;
    pin.(kpv{i}).d=pmd.(kpv{i}).tt;
    pin_out.(kpv{i})=struct2table(pin.(kpv{i}));
    pin_out.(kpv{i})=rmmissing(pin_out.(kpv{i}));
end