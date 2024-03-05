function [plots] = ah_4_edge_plot_6(const, surface, drowh, dpileh, plots)
%I SCREWED UP IN HERE AND EAST/WEST ARE FLOPPED, BUT I JUST CHANGED THE
%FILE NAMES BECAUSE I'M TIRED
dat_bins=const.cott_bins{1};
if const.solve_postgrade==1
    qsurf=surface.Fg;
    bsurf=surface.graded;
    cott=dpileh.cott_pg; % swappable
    pilerev=dpileh.pile_rev_pg; %swappable
else
    qsurf=surface.F_og;
    bsurf=surface.zog;
    cott=dpileh.cotth; % swappable
    pilerev=dpileh.pile_rev; %swappable
end
%create RGB's for dat_bins

yi=200;
xi=13; %should really be an odd number - this used to be 200x21
xpdim=const.xpdim; %panel height in mm
xm(:)=(-xpdim/25.4/2:xpdim/25.4/(xi-1):xpdim/25.4/2)/12;
offset=0.33;%const.offset; %.33 for ati, 2/12 for nevados 1.5/12 NXT
max_tilt=const.maxtilt; %degrees 52 ati 60 nxt
leadedge=const.leading_edge; %leading edge clearance

minx=min(min(surface.xq(~isnan(bsurf))));
maxx=max(max(surface.xq(~isnan(bsurf))));
miny=min(min(surface.yq(~isnan(bsurf))));
maxy=max(max(surface.yq(~isnan(bsurf))));

%% WEST EDGE PLOT
xhc=xm.*sind(max_tilt); %height correction in x
for i=1:height(drowh)
    j=drowh.si(i):drowh.ei(i);
    ph.x{i}=drowh.ntpxc(i)-cosd(max_tilt).*xm; 
    ph.y{i}=drowh.ntpyc(i):-(drowh.rlengthc(i)/yi):drowh.stpyc(i);
    [ph.xm{i}, ph.ym{i}]=meshgrid(ph.x{i},ph.y{i});
    ph.z{i}=interp1(dpileh.tpyc(j),dpileh.tpzc(j),ph.y{i},"linear","extrap");
    ph.zm{i}=ph.z{i}(:)+xhc;
    ph.em{i}=ph.zm{i}-qsurf(ph.xm{i},ph.ym{i});
    ph.gc{i}=ph.em{i}>leadedge;
    ph.rc{i}=ph.em{i}<leadedge;
    ph.bc{i}=ph.rc{i}*0;
end

%PLOT CONTOURS
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 1800])
if (max(max(bsurf)))-min(min(bsurf))<100
    levels=min(min(bsurf))/1000000:1/1000000:max(max(bsurf))/1000000; %1' contours /1000000 so they are under the rows
elseif (max(max(bsurf))-min(min(bsurf)))<200
    levels=min(min(bsurf))/1000000:2/1000000:max(max(bsurf))/1000000;
else
    levels=min(min(bsurf))/1000000:5/1000000:max(max(bsurf))/1000000;
end
cl=contour3(surface.xq,surface.yq,bsurf/1000000,levels,'EdgeColor',[160/255 160/255 160/255]);
hold on

for i=1:numel(ph.x)
    X=reshape(ph.xm{i},[],1);
    Y=reshape(ph.ym{i},[],1);
    Z=reshape(ph.em{i},[],1);
    Z(Z==0)=NaN;
    Z(isnan(Z))=1;
    Zcol=[reshape(ph.rc{i},[],1),reshape(ph.gc{i},[],1),reshape(ph.bc{i},[],1)];
    scatter3(X,Y,Z,5,Zcol);
end
view(2)
axis off
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
daspect([1 1 1])

for i=1:numel(ph.x)
    eapr(i)=sum(sum((ph.em{i}<leadedge)*(xpdim/25.4/(xi-1))*(drowh.rlengthc(i)/yi)));%.*ph.y{i};  %edgegradeareaperrow
    ecpr(i)=sum(sum((ph.em{i}<leadedge).*(ph.em{i}-leadedge)*(xpdim/25.4/(xi-1))*(drowh.rlengthc(i)/yi))); %cut per row
end
formatSpec1 = '%.2f';
formatSpec2 = '%.1f';
wega=num2str(sum(eapr)/43560,formatSpec2); %west edge graded area
wegv=num2str(sum(ecpr)/27,formatSpec1); %west edge graded volume
%dim = [0.3 0.5 0.3 0.3];
%str = {['Total Cut = ' num2str(wegv) ' CY'],['Total Graded Area = ' num2str(wega) ' Acres']};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1,1,0]);
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
saveas(gcf,append(const.fpath{1},'/eastedge.png'))


import mlreportgen.ppt.*
t=string(datetime("today"));

pbt = table(string(wega),string(wegv),'VariableNames',["Graded Area (Acres)","Cut (CY)"]);
ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
    '_EAST_EDGE_INTER_', char(t), '_h.pptx'),'PPT/edge_fig.potx');
eename=(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project), '_EAST_EDGE_INTER_', char(t), '_h.pptx'));
open(ppt);
PMN=Picture(append(const.fpath{1},'/eastedge.png'));
replace(ppt,'PlanarMain',PMN)
%grading brackets
formatSpec = ' PROJECT HAS A MINIMUM EDGE CLEARANCE OF %.1f FEET';
gbstrh=sprintf(formatSpec,const.leading_edge);
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


%% EAST EDGE PLOT

xhc=xm.*sind(-max_tilt); %height correction in x
for i=1:height(drowh)
    j=drowh.si(i):drowh.ei(i);
    ph.x{i}=drowh.ntpxc(i)-cosd(max_tilt).*xm; 
    ph.y{i}=drowh.ntpyc(i):-(drowh.rlengthc(i)/yi):drowh.stpyc(i);
    [ph.xm{i}, ph.ym{i}]=meshgrid(ph.x{i},ph.y{i});
    ph.z{i}=interp1(dpileh.tpyc(j),dpileh.tpzc(j),ph.y{i},"linear","extrap");
    ph.zm{i}=ph.z{i}(:)+xhc;
    ph.em{i}=ph.zm{i}-qsurf(ph.xm{i},ph.ym{i});
    ph.gc{i}=ph.em{i}>leadedge;
    ph.rc{i}=ph.em{i}<leadedge;
    ph.bc{i}=ph.rc{i}*0;
end

%PLOT CONTOURS
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 1800])
if (max(max(bsurf)))-min(min(bsurf))<100
    levels=min(min(bsurf))/1000000:1/1000000:max(max(bsurf))/1000000; %1' contours /1000000 so they are under the rows
elseif (max(max(bsurf))-min(min(bsurf)))<200
    levels=min(min(bsurf))/1000000:2/1000000:max(max(bsurf))/1000000;
else
    levels=min(min(bsurf))/1000000:5/1000000:max(max(bsurf))/1000000;
end
cl=contour3(surface.xq,surface.yq,bsurf/1000000,levels,'EdgeColor',[160/255 160/255 160/255]);
hold on

for i=1:numel(ph.x)
    X=reshape(ph.xm{i},[],1);
    Y=reshape(ph.ym{i},[],1);
    Z=reshape(ph.em{i},[],1);
    Z(Z==0)=NaN;
    Z(isnan(Z))=1;
    Zcol=[reshape(ph.rc{i},[],1),reshape(ph.gc{i},[],1),reshape(ph.bc{i},[],1)];
    scatter3(X,Y,Z,5,Zcol);
end
view(2)
axis off
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
daspect([1 1 1])

for i=1:numel(ph.x)
    eapr(i)=sum(sum((ph.em{i}<leadedge)*(xpdim/25.4/(xi-1))*(drowh.rlengthc(i)/yi)));%.*ph.y{i};  %edgegradeareaperrow
    ecpr(i)=sum(sum((ph.em{i}<leadedge).*(ph.em{i}-leadedge)*(xpdim/25.4/(xi-1))*(drowh.rlengthc(i)/yi))); %cut per row
end
formatSpec1 = '%.2f';
formatSpec2 = '%.1f';
wega=num2str(sum(eapr)/43560,formatSpec2); %west edge graded area
wegv=num2str(sum(ecpr)/27,formatSpec1); %west edge graded volume
%dim = [0.3 0.5 0.3 0.3];
%str = {['Total Cut = ' num2str(wegv) ' CY'],['Total Graded Area = ' num2str(wega) ' Acres']};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1,1,0]);
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
saveas(gcf,append(const.fpath{1},'/westedge.png'))


import mlreportgen.ppt.*
t=string(datetime("today"));

pbt = table(string(wega),string(wegv),'VariableNames',["Graded Area (Acres)","Cut (CY)"]);
ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
    '_WEST_EDGE_INTER_', char(t), '_h.pptx'),'PPT/edge_fig.potx');
wename=(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project), '_WEST_EDGE_INTER_', char(t), '_h.pptx'));
open(ppt);
PMN=Picture(append(const.fpath{1},'/westedge.png'));
replace(ppt,'PlanarMain',PMN)
%grading brackets
formatSpec = ' PROJECT HAS A MINIMUM EDGE CLEARANCE OF %.1f FEET';
gbstrh=sprintf(formatSpec,const.leading_edge);
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

plots.fwestedge=append(const.fpath{1},'/westedge.png');
plots.feastedge=append(const.fpath{1},'/eastedge.png');
if ispc
    pptview(wename,'converttopdf')
    plots.owe=strrep(wename,'pptx','pdf');
    pptview(eename,'converttopdf')
    plots.oee=strrep(eename,'pptx','pdf');
else
    plots.owe=wename;
    plots.oee=eename;
end
end



