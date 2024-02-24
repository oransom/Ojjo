function [ppp_out] = pierplot(kpv, pmd, po, const)
%this PPP includes all mover rows only
%block number
%section
%row number
%truss letter
%tpx
%tpy
%reveal height final in cm
%top pile elevation final in feet
%slope percentage - this request for format came from jesse on january 10
%to orcas email
for i=1:length(kpv)
    ppp.(kpv{i}).block=pmd.(kpv{i}).block';
    ppp.(kpv{i}).section=pmd.(kpv{i}).sect;
    ppp.(kpv{i}).row_number=pmd.(kpv{i}).row;
    ppp.(kpv{i}).truss_letter=pmd.(kpv{i}).tl;
    ppp.(kpv{i}).tpx=po.(kpv{i}).tpx;
    ppp.(kpv{i}).tpy=po.(kpv{i}).tpy;
    ppp.(kpv{i}).reveal_ht_ft=po.(kpv{i}).lp;
    ppp.(kpv{i}).reveal_ht_cm=po.(kpv{i}).lp*30.48;
    ppp.(kpv{i}).top_pile_elev_ft=po.(kpv{i}).tpz;
    ppp.(kpv{i}).nsy=po.(kpv{i}).nsy;
    ppp.(kpv{i}).nsy(isnan(ppp.(kpv{i}).nsy))=0;
    ppp.(kpv{i}).bpy=po.(kpv{i}).bpy;
    ppp.(kpv{i}).bpydiff=ppp.(kpv{i}).bpy-ppp.(kpv{i}).nsy;
    ppp.(kpv{i}).slope_pct=tand(po.(kpv{i}).slp)*100;
    if const.Nevados==1
    ppp.(kpv{i}).NevJoint=po.(kpv{i}).NevJoint;
    ppp.(kpv{i}).LocalSlope=po.(kpv{i}).pierSumSlope;
    end
    ppp_out.(kpv{i})=struct2table(ppp.(kpv{i}));
    ppp_out.(kpv{i})=rmmissing(ppp_out.(kpv{i}));
end
