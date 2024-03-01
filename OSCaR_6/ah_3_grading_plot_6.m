function [plots] = ah_3_grading_plot_6(const, surface, plots) 

minx=min(min(surface.xq(~isnan(surface.zog))));
maxx=max(max(surface.xq(~isnan(surface.zog))));
miny=min(min(surface.yq(~isnan(surface.zog))));
maxy=max(max(surface.yq(~isnan(surface.zog))));

grading.fill=sum(sum(surface.sdiff(surface.sdiff > 0)))*(const.bin_x*const.bin_y)/27; %fill
grading.cut=abs(sum(sum(surface.sdiff(surface.sdiff < 0)))*(const.bin_x*const.bin_y)/27); %cut
grading.area=(sum(sum(surface.sdiff>0))+abs(sum(sum(surface.sdiff<0))))*(const.bin_x*const.bin_y)/43560;

figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 2400])
redChannel = surface.sdiff < 0; %this is cut
greenChannel = surface.sdiff > 0; %this is fill
blueChannel = zeros(size(greenChannel));
colors = double(cat(3, redChannel, greenChannel, blueChannel));
if (max(max(surface.graded(:,:)))-min(min(surface.graded(:,:))))<100
    levels=10000+min(min(surface.graded(:,:)))/1000000:1/1000000:10000+max(max(surface.graded(:,:)))/1000000; %1' contours /1000000 so they are under the rows
elseif (max(max(surface.graded(:,:)))-min(min(surface.graded(:,:))))<200
    levels=10000+min(min(surface.graded(:,:)))/1000000:2/1000000:10000+max(max(surface.graded(:,:)))/1000000;
else
    levels=10000+min(min(surface.graded(:,:)))/1000000:5/1000000:10000+max(max(surface.graded(:,:)))/1000000;
end
contour3(surface.xq,surface.yq,10000+surface.graded/1000000,levels,'-k')
alpha(0.4)
hold on
s=surf(surface.xq, surface.yq, surface.sdiff+10000, colors);
shading interp
view(2)
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
s.AlphaData = (surface.sdiff<=-0.001)|(surface.sdiff>=0.001);
s.FaceAlpha = 0.35;
s.EdgeAlpha = 0.5;
daspect([1 1 1])
axis off
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
saveas(gcf,append(const.fpath{1},'/cutandfill.png'))


formatSpec = '%.1f';
tcut = num2str(grading.cut,formatSpec);
tfill = num2str(grading.fill,formatSpec);
tarea = num2str(grading.area,formatSpec);

import mlreportgen.ppt.*
t=string(datetime("today"));
pbt = table(string(tcut),string(tfill),string(tarea),'VariableNames',["Cut (CY)","Fill (CY)","Disturbed Area (Acres)"]);
if contains(const.tracker,'XTR','IgnoreCase',true)
    ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
        '_GRADING_ESTIMATES_', char(t), '.pptx'),'PPT/grading_figxtr.potx');
else
    ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
        '_GRADING_ESTIMATES_', char(t), '.pptx'),'PPT/grading_fig.potx');
end
gname=(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project), '_GRADING_ESTIMATES_', char(t), '.pptx'));
open(ppt);
PMN=Picture(append(const.fpath{1},'/cutandfill.png'));
replace(ppt,'PlanarMain',PMN)
%grading brackets
formatSpec = ' PROJECT IS GRADED TO MAINTAIN REVEAL RANGE OF %.1f FEET TO %.1f FEET ';
gbstrh=sprintf(formatSpec,const.min_wp,const.max_wp);
gbstr=Paragraph(gbstrh);
gbstr.Style ={HAlign('center')};
gbstr.Font="Hack";
gbstr.FontSize="8";
replace(ppt,"GradeBracket",gbstr)
%project name
formatSpec = 'PROJECT\n%s';
pj=upper(string(const.project));
projnameh=sprintf(formatSpec, pj);
projname = Paragraph(projnameh);
projname.Style = {HAlign('center')};
projname.Font="Hack";
projname.FontSize="12";
replace(ppt,"ProjectCoSt",projname)
%customer
formatSpec = '%s';
cn=upper(string(const.customer));
custnameh=sprintf(formatSpec, cn);
custname = Paragraph(custnameh);
custname.Style = {HAlign('center')};
custname.Font="Hack";
custname.FontSize="12";
replace(ppt,"Customer",custname)
%customer
formatSpec = '%s\n%s';
city=upper(string(const.city));
co_st=upper(string(const.county));
citycoh=sprintf(formatSpec, city, co_st);
cityco = Paragraph(citycoh);
cityco.Style = {HAlign('center')};
cityco.Font="Hack";
cityco.FontSize="12";
replace(ppt,"Coordinates",cityco)


%revision table
revtable=Table({'Rev','Date','By','Checked','Approved','Description';...
    const.rev_num, t, "OR", "JR", "OR","AutoGrade";...
    " ", " ", " ", " "," "," "});
revtable.Font="Hack";
revtable.FontSize="6";
revtable.entry(1,1).Style={ColWidth('.35in')};
revtable.entry(1,2).Style={ColWidth('.7in')};
revtable.entry(1,3).Style={ColWidth('.30in')};
revtable.entry(1,4).Style={ColWidth('.5in')};
revtable.entry(1,5).Style={ColWidth('.55in')};
revtable.entry(1,6).Style={ColWidth('1.02in')};
revtable.row(1).Style = [revtable.row(1).Style {RowHeight("0.06in")}];
revtable.row(2).Height="0.06in";
revtable.row(3).Height="0.06in";
replace(ppt,'RevTable',revtable)

%PileTable
PileTable=Table(pbt);
PileTable.Font="Hack";
PileTable.FontSize="6";
PileTable.entry(1,1).Style={ColWidth('0.8in')};
PileTable.entry(1,2).Style={ColWidth('0.8in')};
PileTable.entry(1,2).Style={ColWidth('0.8in')};
PileTable.entry(1,2).Style={ColWidth('0.8in')};
PileTable.row(1).Style = [PileTable.row(1).Style {RowHeight("0.06in")}];
for j=2:size(pbt,1)
    PileTable.row(j).Height="0.06in";
end
replace(ppt,'PileTable',PileTable)
close(ppt);

plots.fgrding=append(const.fpath{1},'/cutandfill.png');
if ispc
    pptview(gname,'converttopdf')
    plots.grading=strrep(gname,'pptx','pdf');
else
    plots.grading=gname;
end

end