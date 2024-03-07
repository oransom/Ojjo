function [plots] = ag_1_dataRollup_6(const,drowh,dpileh,grading,plots)

  

writetable(dpileh,append(const.dpath{1},'/Pile_Data_All.xlsx'))
writetable(drowh,append(const.dpath{1},'/Row_Data_All.xlsx'))

ru.cut=round(grading.cut,0);
ru.fill=round(grading.fill,0);
ru.grdara=round(grading.area,2);

ru.pregrdagstl=round(sum(dpileh.pile_rev),1);
if ~contains(const.tracker,'ojjo','IgnoreCase',true)
    ru.pregrddrvstl=round(sum(dpileh.driven_steel),1);
    ru.pregrdtotstl=round(sum(dpileh.total_steel),1);


    if const.solve_postgrade==1

        ru.postgrdagstl=round(sum(dpileh.pile_rev_pg),1);
        ru.postgrddrvstl=round(sum(dpileh.driven_steel_pg),1);
        ru.postgrdtotstl=round(sum(dpileh.total_steel_pg),1);
    end
end

ru=struct2table(ru);
writetable(ru,append(const.rupath{1},'/',const.project{1},'_Rollup.xlsx'))
end