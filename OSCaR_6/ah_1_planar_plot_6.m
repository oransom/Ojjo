function [plots] = ah_1_planar_plot_6(const, surface, drowh, dpileh)
%change some variable names for compactness - not best for memory, but
%whatever
%tic
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
for j=1:length(dat_bins.bin_color)
    if matches(string(dat_bins.bin_color(j)),"b")
        dat_bins.rgb{j}=[0 0 1];
    elseif matches(string(dat_bins.bin_color(j)),"lb")
        dat_bins.rgb{j}=[.8 1 1];
    elseif matches(string(dat_bins.bin_color(j)),"db")
        dat_bins.rgb{j}=[0 0 .5];
    elseif matches(string(dat_bins.bin_color(j)),"g")
        dat_bins.rgb{j}=[0 1 0];
    elseif matches(string(dat_bins.bin_color(j)),"lg")
        dat_bins.rgb{j}=[.78 1 .78];
    elseif matches(string(dat_bins.bin_color(j)),"dg")
        dat_bins.rgb{j}=[0 .42 0];
    elseif matches(string(dat_bins.bin_color(j)),"y")
        dat_bins.rgb{j}=[1 1 0];
    elseif matches(string(dat_bins.bin_color(j)),"o")
        dat_bins.rgb{j}=[1 .4 0];
    elseif matches(string(dat_bins.bin_color(j)),"do")
        dat_bins.rgb{j}=[.42 .21 0];
    elseif matches(string(dat_bins.bin_color(j)),"lo")
        dat_bins.rgb{j}=[1 .898 .8];
    elseif matches(string(dat_bins.bin_color(j)),"r")
        dat_bins.rgb{j}=[1 0 0];
    elseif matches(string(dat_bins.bin_color(j)),"lr")
        dat_bins.rgb{j}=[1 .6 .6];
    elseif matches(string(dat_bins.bin_color(j)),"dr")
        dat_bins.rgb{j}=[.6 0 0];
    else
        dat_bins.rgb{j}=[1 0 0];
    end
end

yi=200;
xi=13; %should really be an odd number - this used to be 200x21
xpdim=const.xpdim; %panel height in mm
xm(:)=(-xpdim/25.4/2:xpdim/25.4/(xi-1):xpdim/25.4/2)/12;
offset=0;%const.offset; %.33 for ati, 2/12 for nevados 1.5/12 NXT
max_tilt=const.maxtilt; %degrees 52 ati 60 nxt
leadedge=const.leading_edge; %leading edge clearance

minx=min(min(surface.xq(~isnan(bsurf))));
maxx=max(max(surface.xq(~isnan(bsurf))));
miny=min(min(surface.yq(~isnan(bsurf))));
maxy=max(max(surface.yq(~isnan(bsurf))));

for i=1:height(drowh)
    j=drowh.si(i):drowh.ei(i);
    ph.x{i}=drowh.ntpxc(i)-xm;
    ph.y{i}=drowh.ntpyc(i):-(drowh.rlengthc(i)/yi):drowh.stpyc(i);
    [ph.xm{i}, ph.ym{i}]=meshgrid(ph.x{i},ph.y{i});
    ph.z{i}=interp1(dpileh.tpyc(j),dpileh.tpzc(j),ph.y{i},"linear","extrap");
    ph.zm{i}=repmat(ph.z{i}',[1,xi]);
    ph.em{i}=ph.zm{i}-qsurf(ph.xm{i},ph.ym{i});
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
%PLANAR PLOTTER
for k=1:length(dat_bins.bin_start)
    if k==1
        binz(k,1:2)=[-1000,dat_bins.bin_finish(k)];
    elseif k==length(dat_bins.bin_start)
        binz(k,1:2)=[dat_bins.bin_start(k),1000];
    else
        binz(k,1:2)=[dat_bins.bin_start(k),dat_bins.bin_finish(k)];
    end
end
binzhc=reshape(binz,1,[]);
binzhc=sort(binzhc);
binzhc=unique(binzhc);
for i=1:numel(ph.x)
    X=reshape(ph.xm{i},[],1);
    Y=reshape(ph.ym{i},[],1);
    Z=reshape(ph.em{i},[],1);
    Z(Z==0)=NaN;
    Z(isnan(Z))=1;
    Zc=discretize(Z,binzhc);
    Zcol=cell2mat(dat_bins.rgb(Zc));
    scatter3(X,Y,Z,5,Zcol);
end
daspect([1 1 1])
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [    for j=1:length(fliped.tpxc)X Y W H]
set(gca, 'Position', NewPos);

%%GRADING EXTENT PLOTTER
% if ~isempty(surface.gbound)
%     for j=1:numel(surface.gbound)
%         plot3(surface.xg{j}(surface.gbound{j}),...
%               surface.yg{j}(surface.gbound{j}),...
%               (ones(numel(surface.gbound{j}),1).*1000000),'k-','LineWidth',2)
%         hold on
%     end
% end
if ~isempty(surface.gbound)
    for j=1:numel(surface.gbound)
        hx=patch('XData',surface.xg{j}(surface.gbound{j}),...
                 'YData',surface.yg{j}(surface.gbound{j}),...
                 'ZData',(ones(numel(surface.gbound{j}),1)*(max(cott)+10)),...
                 'EdgeColor','k','LineWidth',1);
        hatchfill2(hx,'cross','LineWidth',1,'FaceColor','none','HatchStyle','single','HatchAngle',135,'HatchDensity',300);
    end
end
view(2)
axis off
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])

%FLIPEX
if strcmpi('flip2ext',drowh.Properties.VariableNames)
box_x=xpdim/304.8+.5;
flipex.tpxc=drowh.mdptx(drowh.flip2ext==1);
flipex.tpyc=drowh.mdpty(drowh.flip2ext==1);
flipex.tpzc=drowh.rowzavg(drowh.flip2ext==1);
flipex.row_len=drowh.rlengthc(drowh.flip2ext==1);

if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
    fliped.tpxc=drowh.mdptx(drowh.flip2edg==1);
    fliped.tpyc=drowh.mdpty(drowh.flip2edg==1);
    fliped.tpzc=drowh.rowzavg(drowh.flip2edg==1);
    fliped.row_len=drowh.rlengthc(drowh.flip2edg==1);
    for j=1:length(fliped.tpxc)
        box_y=fliped.row_len(j)/2+2;
        P = [fliped.tpxc(j)-box_x fliped.tpyc(j)-box_y fliped.tpzc(j); fliped.tpxc(j)-box_x fliped.tpyc(j)+box_y fliped.tpzc(j);...
            fliped.tpxc(j)+box_x fliped.tpyc(j)+box_y fliped.tpzc(j); fliped.tpxc(j)+box_x fliped.tpyc(j)-box_y fliped.tpzc(j);...
            fliped.tpxc(j)-box_x fliped.tpyc(j)-box_y fliped.tpzc(j)];
        plot3(P(:,1),P(:,2),P(:,3),'r--', 'linewidth',1.5)
    end
    for j=1:length(flipex.tpxc)
        box_y=flipex.row_len(j)/2+2;
        P = [flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j); flipex.tpxc(j)-box_x flipex.tpyc(j)+box_y flipex.tpzc(j);...
            flipex.tpxc(j)+box_x flipex.tpyc(j)+box_y flipex.tpzc(j); flipex.tpxc(j)+box_x flipex.tpyc(j)-box_y flipex.tpzc(j);...
            flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j)];
        plot3(P(:,1),P(:,2),P(:,3),'r-', 'linewidth',1.5)
    end
else
    for j=1:length(flipex.tpxc)
        box_y=flipex.row_len(j)/2+2;
        P = [flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j); flipex.tpxc(j)-box_x flipex.tpyc(j)+box_y flipex.tpzc(j);...
            flipex.tpxc(j)+box_x flipex.tpyc(j)+box_y flipex.tpzc(j); flipex.tpxc(j)+box_x flipex.tpyc(j)-box_y flipex.tpzc(j);...
            flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j)];
        plot3(P(:,1),P(:,2),P(:,3),'r-', 'linewidth',1.5)
    end
end
end

%% FLOODPLOT
if strcmpi(const.flood,'na')==0
    fhold=surface.flood;
    fhold(surface.floodmask)=NaN;
    fhold(isnan(bsurf))=NaN;
    fhold=fhold-bsurf;
    fhold(fhold<0.1)=0;
    fs=surf(surface.xq,surface.yq,fhold);
    clim([min(min(fhold)),max(max(fhold))])
    shading interp
    alpha(fs,0.5);
    c=colorbar("Location","southoutside");
    c.Label.String = 'Flood depth in feet';
end

view(2)
axis off
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
saveas(gcf,append(const.fpath{1},'/heightplot.png'))

%CREATE TABLE NAMES
for k=1:height(dat_bins)
    if k==1
        epname{k}=[num2str(dat_bins.bin_start(k)+const.mtr_cott2top,'%.2f'),'ft - ',num2str(dat_bins.bin_finish(k)+const.std_cott2top,'%.2f'),'ft'];
    elseif k<height(dat_bins)
        epname{k}=[num2str(dat_bins.bin_start(k)+const.std_cott2top,'%.2f'),'ft - ',num2str(dat_bins.bin_finish(k)+const.std_cott2top,'%.2f'),'ft'];
    elseif k==height(dat_bins)
        if dat_bins.bin_start(k)+const.std_cott2top<max(cott)
            epname{k}=[num2str(dat_bins.bin_start(k)+const.std_cott2top,'%.2f'),'ft - ',num2str(max(cott),'%.2f'),'ft'];
        else
            epname{k}=[num2str(dat_bins.bin_start(k)+const.std_cott2top,'%.2f'),'ft +'];
        end
    end
end
e1name=string(epname'); %this is ATI names
edges_pile=[vertcat(dat_bins.bin_start(1)+const.mtr_cott2top,dat_bins.bin_finish(1:end)+const.std_cott2top)]';
edges2=[vertcat(dat_bins.bin_start(1),dat_bins.bin_finish(1:end))]';
e2name=[string(dat_bins.bin_name)]; %this is NXT names

%STANDALONE COLOR BAR
fig_leg=figure('Visible', 'off'); %create stand alone color bar
axis off
for k=1:length(dat_bins.bin_start)
    fig_leg=scatter(nan,nan, 'MarkerFaceColor', dat_bins.rgb{k});
    hold on
end
legend(dat_bins.bin_name')
% Initial values to capture the entire legend
% Should fit most modern screens
set(gcf,'Position',[0,0,100,100]);
% Call the legend to your choice, I used a horizontal legend here
legend_handle = legend('Orientation','vertical');
% Set the figure Position using the normalized legend Position vector
% as a multiplier to the figure's current position in pixels
% This sets the figure to have the same size as the legend
set(gcf,'Position',(get(legend_handle,'Position')...
    .*[0, 0, 1, 1].*get(gcf,'Position')));
% The legend is still offset so set its normalized position vector to
% fill the figure
set(legend_handle,'Position',[0,0,1,1]);
% Put the figure back in the middle screen area
set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);
saveas(gcf,append(const.fpath{1},'/heightcolorbar.png'))

%HEIGHT ABOVE GROUND HISTOGRAM
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 600, 400])
edges=dat_bins.bin_start(1):.25:dat_bins.bin_finish(end);

if max(cott)+.01>edges2(end)
    edges2(end)=max(cott)+.01;
end
if max(cott)>edges_pile(end)
    edges_pile(end)=max(cott);
end
%% PBINS AND CBINS
[pbin, ~]=histcounts(pilerev,edges_pile+.1);
[cottbin, ~]=histcounts(cott,edges2+.1);
pbinp=100*(pbin/(sum(pbin(:))));
cbinp=100*(cottbin/(sum(cottbin(:))));
formatSpec = '%.1f';
for j=1:size(pbinp,2)
    pbins{j} = num2str(pbinp(j),formatSpec);
    cbins{j} = num2str(cbinp(j),formatSpec);
end
if strcmpi(string(const.tracker),"NXT")==1 || contains(string(const.tracker),"XTR",'IgnoreCase',true)
    pbt = table(e1name,pbins','VariableNames',["Pile Reveal","Percentage"]);
else
    pbt = table(e1name,pbins',e2name,cbins','VariableNames',["Top Pile Height","Percentage","CoTT Height","Percentage "]);
end
%%
h1=histogram(cott,edges,'Normalization','probability');
yl=(get(gca,'YLim'));
set(gca,'YLim',[yl(1),yl(2)+.13*yl(2)]);
ylabel('Percentage of Panels')
xlabel('Panel Height Above Ground')
n=get(h1,'Values');
xn=get(h1,"BinEdges");
formatSpec = '%.2f';
barstrings=num2str(n',formatSpec);
text(xn(2:end), n, barstrings,'horizontalalignment','center','verticalalignment','bottom','Rotation',70)
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
saveas(gcf,append(const.fpath{1},'/height_hist.png'))


%% CREATE POWERPOINT SLIDE
import mlreportgen.ppt.*
t=string(datetime("today"));
% Create Slides
switch char(const.tracker)
    case {'ATI'}
        ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_' ,char(const.project),...
            '_Eng_RPCS_Planar_', char(t), '.pptx'),'PPT/RPCS Planar Template - ATI.potx');
        fname=append(const.fpath{1}, '/', char(const.customer), '_' ,char(const.project),...
            '_Eng_RPCS_Planar_', char(t), '.pptx');
    case {'NXT','XTR','XTR1p5'}
        ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
            '_Eng_RPCS_Planar_', char(t), '.pptx'),'PPT/RPCS Planar Template - NXT.potx');
        fname=append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
            '_Eng_RPCS_Planar_', char(t), '.pptx');
    case {'NEV','Ojjo_ATI','Ojjo_NXT','FTC'}
        ppt = Presentation(append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
            '_Eng_Study_', char(t), '.pptx'),'PPT/PPT Template - NEV.potx');
        fname=append(const.fpath{1}, '/', char(const.customer), '_', char(const.project),...
            '_Eng_Study_', char(t), '.pptx');        
end
open(ppt);

PMN=Picture(append(const.fpath{1},'/heightplot.png'));
HST=Picture(append(const.fpath{1},'/height_hist.png'));
CLB=Picture(append(const.fpath{1},'/heightcolorbar.png'));
FLP=Picture('PPT/flipex.png');
if strcmpi('flip2ext',drowh.Properties.VariableNames)
    if sum(drowh.flip2ext)>0
        replace(ppt,'Flipex',FLP)
    end
    FLPD=Picture('PPT/fliped.png');
    if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
        if sum(drowh.flip2edg)>0
            replace(ppt,'Fliped',FLPD)
        end
    end
end
GRD=Picture('PPT/graded.png');
if ~isempty(surface.gbound)
    replace(ppt,'Grading',GRD)
end
replace(ppt,'PlanarMain',PMN)
replace(ppt,'Histogram',HST)
replace(ppt,'ColorBar',CLB)
%max min slope
formatSpec = ' APPROXIMATE MAX NORTH ROW SLOPE = %.1f%c \n APPROXIMATE MAX SOUTH ROW SLOPE = %.1f%c ';
northslope=-min(drowh.slp(drowh.slp<0));
if isempty(northslope)==1
    northslope=0;
end
southslope=max(drowh.slp(drowh.slp>0));
if isempty(southslope)==1
    southslope=0;
end
slpstrh=sprintf(formatSpec, northslope, char(176), southslope, char(176));
slpstr = Paragraph(slpstrh);
slpstr.Font="Hack";
slpstr.FontSize="8";
replace(ppt,"RowSlope",slpstr)
%max min tth
if strcmpi(string(const.tracker),"NXT")==1 || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
    formatSpec = ' MAX PIER HEIGHT = %.1f ft \n MIN PIER HEIGHT = %.1f ft';
    tthstrh=sprintf(formatSpec, max(cott), min(cott));
else
    formatSpec = ' MAX TORQUE TUBE HEIGHT = %.1f ft \n MIN TORQUE TUBE HEIGHT = %.1f ft';
    tthstrh=sprintf(formatSpec, max(cott), min(cott));
end
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
%topo source
formatSpec = 'TOPO SOURCE:\n%s';
ts=upper(string(const.topo_source));
toponameh=sprintf(formatSpec, ts);
toponame = Paragraph(toponameh);
toponame.Style = {HAlign('left')};
toponame.Font="Hack";
toponame.FontSize="12";
replace(ppt,"TopoSource",toponame)
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
    const.rev_num, t, "OR", "JR", "HC","AutoPlanar";...
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
if strcmpi(string(const.tracker),"NXT")==1  || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
    PileTable=Table(pbt);
    PileTable.Font="Hack";
    PileTable.FontSize="6";
    PileTable.entry(1,1).Style={ColWidth('1.6in')};
    PileTable.entry(1,2).Style={ColWidth('1.6in')};
    PileTable.row(1).Style = [PileTable.row(1).Style {RowHeight("0.06in")}];
    for j=2:size(pbt,1)
        PileTable.row(j).Height="0.06in";
    end
    replace(ppt,'PileTable',PileTable)
else
    PileTable=Table(pbt);
    PileTable.Font="Hack";
    PileTable.FontSize="6";
    PileTable.entry(1,1).Style={ColWidth('1.0in')};
    PileTable.entry(1,2).Style={ColWidth('0.6in')};
    PileTable.entry(1,2).Style={ColWidth('1.0in')};
    PileTable.entry(1,2).Style={ColWidth('0.6in')};
    PileTable.row(1).Style = [PileTable.row(1).Style {RowHeight("0.06in")}];
    for j=2:size(pbt,1)
        PileTable.row(j).Height="0.06in";
    end
    replace(ppt,'PileTable',PileTable)
end
close(ppt);

plots.fhghtplot=append(const.fpath{1},'/heightplot.png');
plots.fhghthist=append(const.fpath{1},'/height_hist.png');
plots.fheightcb=append(const.fpath{1},'/heightcolorbar.png');
if ispc
    pptview(planar.fname,'converttopdf')
    plots.planar=strrep(fname,'pptx','pdf');
else
    plots.planar=fname;
end
end



