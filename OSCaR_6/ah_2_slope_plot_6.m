function [plots] = ah_2_slope_plot_6(const, surface, drowh, dpileh, plots)
%change some variable names for compactness - not best for memory, but
%whatever

qsurf=surface.F_og;
bsurf=surface.zog;
cott=dpileh.cotth; % swappable
pilerev=dpileh.pile_rev; %swappable
box_x=const.xpdim/304.8+.5;
minx=min(min(surface.xq(~isnan(bsurf))));
maxx=max(max(surface.xq(~isnan(bsurf))));
miny=min(min(surface.yq(~isnan(bsurf))));
maxy=max(max(surface.yq(~isnan(bsurf))));

%PLOT CONTOURS
figure ('Visible', 'off')
hold on
set(gcf, 'Position',  [100, 100, 2400, 1800])
if (max(max(bsurf)))-min(min(bsurf))<100
    levels=min(min(bsurf))/1000000:1/1000000:max(max(bsurf))/1000000; %1' contours /1000000 so they are under the rows
elseif (max(max(bsurf))-min(min(bsurf)))<200
    levels=min(min(bsurf))/1000000:2/1000000:max(max(bsurf))/1000000;
else
    levels=min(min(bsurf))/1000000:5/1000000:max(max(bsurf))/1000000;
end
cl=contour3(surface.xq,surface.yq,bsurf/1000000,levels,'EdgeColor',[160/255 160/255 160/255]);

%SLOPE BOXES
for i=1:height(drowh)
    if max(abs(drowh.slp(i)))<4.0
        Pz4=[min(drowh.mdptx(i))-box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))-box_x max(drowh.mdpty(i))+(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))+box_x max(drowh.mdpty(i))+(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))+box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))-box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i))];
        patch('XData',Pz4(:,1),'YData',Pz4(:,2),'ZData',Pz4(:,3),'EdgeColor','g','FaceColor','g',...
              'FaceAlpha',0.25,'EdgeAlpha',0.5)
    elseif max(abs(drowh.slp(i)))>4.0 && max(abs(drowh.slp(i)))<8.0
        P48=[min(drowh.mdptx(i))-box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))-box_x max(drowh.mdpty(i))+(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))+box_x max(drowh.mdpty(i))+(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))+box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))-box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i))];
        patch('XData',P48(:,1),'YData',P48(:,2),'ZData',P48(:,3),'EdgeColor', [1 .6412 0],'FaceColor', [1 .6412 0],...
              'FaceAlpha',0.25,'EdgeAlpha',0.5)
    else
        P8P=[min(drowh.mdptx(i))-box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))-box_x max(drowh.mdpty(i))+(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))+box_x max(drowh.mdpty(i))+(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))+box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2) max(drowh.rowzavg(i));...
            min(drowh.mdptx(i))-box_x min(drowh.mdpty(i))-(drowh.rlengthc(i)/2+3) max(drowh.rowzavg(i))];
        patch('XData',P8P(:,1),'YData',P8P(:,2),'ZData',P8P(:,3),'EdgeColor', 'r','FaceColor', 'r',...
              'FaceAlpha',0.25,'EdgeAlpha',0.5)
    end
end


view(2)
axis off
hold on
xlim([min(minx)-10 max(maxx)+10])
ylim([min(miny)-10 max(maxy)+10])

daspect([1 1 1])
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
saveas(gcf,append(const.fpath{1},'/rowslope_plot.png'))

formatSpec = '%.1f';
slpm=abs(drowh.slp);
sl4=100*(sum(slpm<4.0)/length(slpm));
s4t8=100*(sum(slpm>4.0 & slpm<8.0)/length(slpm));
sg8=100*(sum(slpm>8.0)/length(slpm));
tcut = num2str(sl4,formatSpec);
tfill = num2str(s4t8,formatSpec);
tarea = num2str(sg8,formatSpec);

import mlreportgen.ppt.*
t=string(datetime("today"));
pbt = table(string(tcut),string(tfill),string(tarea),'VariableNames',["Rows 0˚-4.0˚","Rows 4.0˚-8.0˚","Rows > 8.0˚"]);

ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_' ,char(const.project),...
    '_ROW_SLOPE_', char(t), '.pptx'),'PPT/rowslope_fig.potx');
sname=append(const.fpath{1}, '/', char(const.customer), '_' ,char(const.project),...
    '_ROW_SLOPE_', char(t), '.pptx');
open(ppt);
PMN=Picture(append(const.fpath{1},'/rowslope_plot.png'));
replace(ppt,'PlanarMain',PMN)

if min(drowh.slp)<0
    nslp=-min(drowh.slp);
else
    nslp=0;
end
if max(drowh.slp>0)
    sslp=max(drowh.slp);
else
    sslp=0;
end


%max min slope
formatSpec = ' APPROXIMATE MAX NORTH ROW SLOPE = %.1f%c \n APPROXIMATE MAX SOUTH ROW SLOPE = %.1f%c ';
slpstrh=sprintf(formatSpec, nslp, char(176), sslp, char(176));
slpstr = Paragraph(slpstrh);
slpstr.Font="Hack";
slpstr.FontSize="8";
replace(ppt,"RowSlope",slpstr)
%max min tth
formatSpec = ' MAX TORQUE TUBE HEIGHT = %.1f ft \n MIN TORQUE TUBE HEIGHT = %.1f ft';
tthstrh=sprintf(formatSpec, max(cott), min(cott));
tthstr = Paragraph(tthstrh);
tthstr.Font="Hack";
tthstr.FontSize="8";
replace(ppt,"MinMaxRow",tthstr)
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
    const.rev_num, t, "OR", "JR", "OR","AutoRowSlope";...
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
PileTable.entry(1,1).Style={ColWidth('1.0in')};
PileTable.entry(1,2).Style={ColWidth('1.0in')};
PileTable.entry(1,2).Style={ColWidth('1.0in')};
%PileTable.entry(1,2).Style={ColWidth('0.8in')};
PileTable.row(1).Style = [PileTable.row(1).Style {RowHeight("0.06in")}];
for j=2:size(pbt,1)
    PileTable.row(j).Height="0.06in";
end
replace(ppt,'PileTable',PileTable)
close(ppt);
plots.fslpplot=append(const.fpath{1},'/rowslope_plot.png');
if ispc
    pptview(sname,'converttopdf')
    plots.slope=strrep(sname,'pptx','pdf');
else
    plots.slope=sname;
end


end