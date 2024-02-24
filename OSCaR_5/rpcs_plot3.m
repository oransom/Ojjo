function [fname, fnamepdf, pbt] = rpcs_plot3(kpv,pd_out,F_og,xyzf_int,surface,ppp_out,const,grd_ext,fliped,flipex)
%change some variable names for compactness - not best for memory, but
%whatever
%tic 
p=ppp_out; %tpx=zeros(200,200);
dat_bins=const.cott_bins{1};
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

%rollup some data
for i=1:length(kpv)
    if i==1
        cott_ft_ru=pd_out.(kpv{i}).cott_ft;
        pile_reveal_ru=pd_out.(kpv{i}).pile_reveal;
        slp_deg_ru=pd_out.(kpv{i}).slp_deg;
    else
        cott_ft_ru=[cott_ft_ru;pd_out.(kpv{i}).cott_ft];
        pile_reveal_ru=[pile_reveal_ru;pd_out.(kpv{i}).pile_reveal];
        slp_deg_ru=[slp_deg_ru;pd_out.(kpv{i}).slp_deg];
    end
end

cott_ft_ru=floor(cott_ft_ru*10)/10;
pile_reveal_ru=floor(pile_reveal_ru*10)/10;


figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 1800])
if (max(max(surface.zog)))-min(min(surface.zog))<100
    levels=min(min(surface.zog))/1000000:1/1000000:max(max(surface.zog))/1000000; %1' contours /1000000 so they are under the rows
elseif (max(max(surface.zog))-min(min(surface.zog)))<200
    levels=min(min(surface.zog))/1000000:2/1000000:max(max(surface.zog))/1000000;
else
    levels=min(min(surface.zog))/1000000:5/1000000:max(max(surface.zog))/1000000;
end
contour3(surface.xq,surface.yq,surface.zog/1000000,levels,'-k')
alpha(0.4)
hold on
%% ANGLE == FLAT
for i=1:length(kpv)
    numsect=1:max(p.(kpv{i}).section);
    for j=1:max(numsect)
        if ismember(j,p.(kpv{i}).section)
            numrow(j)=max(p.(kpv{i}).row_number(p.(kpv{i}).section==j));
        else
            numrow(j)=NaN;
        end
    end
    l=0;
    for j=1:max(numsect)
        for k=1:max(numrow(j))
            l=l+1;
            tpx{:,l}=p.(kpv{i}).tpx(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
            tpy{:,l}=p.(kpv{i}).tpy(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
            top{:,l}=p.(kpv{i}).reveal_ht_ft(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
            elev{:,l}=p.(kpv{i}).top_pile_elev_ft(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
        end
    end
    minx(i)=min(p.(kpv{i}).tpx(:));
    maxx(i)=max(p.(kpv{i}).tpx(:));
    miny(i)=min(p.(kpv{i}).tpy(:));
    maxy(i)=max(p.(kpv{i}).tpy(:));

    %if a error happens somewhere around here check that each section is restarting at 1 and not continuing
    for j=1:length(tpx) %parfor
        if ~isempty(tpx{j})
            x=cell2mat(tpx(j));
            y=cell2mat(tpy(j));
            h=cell2mat(top(j));
            e=cell2mat(elev(j));
            y_int=repmat(y(1):(y(end)-y(1))/yi:y(end),xi,1);
            yh(:,:,j)=y_int;
            x_int=ones(size(y_int)).*(xm+mean(x))';
            xh(:,:,j)=x_int;
            e_int=interp1(y,e,y_int);
            eh(:,:,j)=e_int;
        end
    end
    for j=1:length(tpx)
        if ~isempty(tpx{j})
            %pre or post grade
            if const.solve_postgrade==1
                hag=eh(:,:,j)-xyzf_int(xh(:,:,j),yh(:,:,j));
            else
                hag=eh(:,:,j)-F_og(xh(:,:,j),yh(:,:,j)); %for future FYI we can use histcounts here with edges to bin distance above ground
            end
            panels.(kpv{i}).hagh(:,:,j)=hag+offset;
        end
    end
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
    for j=1:length(tpx)
        %if ~isempty(tpx{j})
            if j==1 && ~isempty(tpx{j})
                X=reshape(xh(:,:,j),[],1);
                Y=reshape(yh(:,:,j),[],1);
                Z=reshape(panels.(kpv{i}).hagh(:,:,j),[],1);
            elseif j==1 && isempty(tpx{j})
                X=[];
                Y=[];
                Z=[];
            elseif ~isempty(tpx{j})
                X=[X;reshape(xh(:,:,j),[],1)];
                Y=[Y;reshape(yh(:,:,j),[],1)];
                Z=[Z;reshape(panels.(kpv{i}).hagh(:,:,j),[],1)];
            end
        %end
    end
    Z(Z==0)=NaN;
    Z(isnan(Z))=1;
    Zc=discretize(Z,binzhc);
    Zcol=cell2mat(dat_bins.rgb(Zc));
    scatter3(X,Y,Z,5,Zcol);
    diz2=toc;
    fprintf('planar scattering complete %3.2f s \n', diz2);

    
    %GRADING EXTENT PLOTTER
    if ~isempty(grd_ext)
        for j=1:length(grd_ext.gbx)
            plot(grd_ext.xgh{j}(grd_ext.gbx{j}),grd_ext.ygh{j}(grd_ext.gby{j}),'k--','LineWidth',2)
            hold on
        end
    end

    view(2)
    axis off
    hold on
    xlim([min(minx)-50 max(maxx)+50])
    ylim([min(miny)-50 max(maxy)+50])

    daspect([1 1 1])
    set(gca, 'units', 'normalized'); %Just making sure it's normalized
    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
    %[Left Bottom Right Top] spacing
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [    for j=1:length(fliped.tpxc)X Y W H]
    set(gca, 'Position', NewPos);

    if const.min_wp ~= const.max_wp
        clim([const.min_wp const.max_wp])
    else
        clim([const.min_wp+offset-.5 const.max_wp+offset+.5])
    end
end
%flipex boxes
box_x=xpdim/304.8+.5;

if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
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
if strcmpi(const.flood,'na')==0
    maxx=max(max(surface.xq));
    minx=min(min(surface.xq));
    maxy=max(max(surface.yq));
    miny=min(min(surface.yq));
    surface.floodplot=surface.flood;
    %surface.floodplot(surface.floodmask)=NaN;
    surface.floodplot=surface.floodplot-surface.untrimmed;
    fs=surf(surface.xq,surface.yq,surface.floodplot);
    clim([min(min(surface.floodplot)),max(max(surface.floodplot))])
    shading interp
    alpha(fs,0.5);
    c=colorbar("Location","southoutside");
    c.Label.String = 'Flood depth in feet';
end
if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' height_plot.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' height_plot.png'])
    end
    foo=1;
end
%create some table names
for k=1:height(dat_bins)
    if k==1
        epname{k}=[num2str(dat_bins.bin_start(k)-const.mtr_wp2top,'%.2f'),'ft - ',num2str(dat_bins.bin_finish(k)-const.std_wp2top,'%.2f'),'ft'];
    elseif k<height(dat_bins)
        epname{k}=[num2str(dat_bins.bin_start(k)-const.std_wp2top,'%.2f'),'ft - ',num2str(dat_bins.bin_finish(k)-const.std_wp2top,'%.2f'),'ft'];
    elseif k==height(dat_bins)
        epname{k}=[num2str(dat_bins.bin_start(k)-const.std_wp2top,'%.2f'),'ft - ',num2str(max(cott_ft_ru),'%.2f'),'ft'];%pd_out.(kpv{i}).cott_ft
    end
end
e1name=string(epname');

%% Other Plots
%Colorbar for height above ground plot
fig_leg=figure('Visible', 'off'); %create stand alone color bar
axis off
for k=1:length(dat_bins.bin_start)
    fig_leg=scatter(nan,nan, 'MarkerFaceColor', dat_bins.rgb{k});
    hold on
end
if strcmpi(string(const.tracker),"NXT")==1 || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
    legend(dat_bins.bin_name')
else
    legend(dat_bins.bin_name')
end
% Initial values to capture the entire legend
% Should fit most modern screens
set(gcf,'Position',[0,0,1024,1024]);
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

if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png'])
    end
    foo=1;
end

%height above ground histogram
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 600, 400])
edges=[dat_bins.bin_start(1):.25:dat_bins.bin_finish(end)];
%edges=[3.0:.25:10];
edges_pile=[vertcat(dat_bins.bin_start(1)-const.mtr_wp2top,dat_bins.bin_finish(1:end)-const.std_wp2top)]';
edges2=[vertcat(dat_bins.bin_start(1),dat_bins.bin_finish(1:end))]';
e2name=[string(dat_bins.bin_name)];

if max(cott_ft_ru)+.1>edges2(end)
    edges2(end)=max(cott_ft_ru)+.1;
end
if max(cott_ft_ru)>edges_pile(end)
    edges_pile(end)=max(cott_ft_ru);
end
%% pbins and cbins
[pbin pbedges]=histcounts(pile_reveal_ru,edges_pile+.1);%pd_out.(kpv{i}).pile_reveal
[cottbin cottedges]=histcounts(cott_ft_ru,edges2+.1);%pd_out.(kpv{i}).cott_ft
pbinp=100*(pbin/(sum(pbin(:))));
cbinp=100*(cottbin/(sum(cottbin(:))));
formatSpec = '%.1f';
for j=1:size(pbinp,2)
    pbins{j} = num2str(pbinp(j),formatSpec);
    cbins{j} = num2str(cbinp(j),formatSpec);
end
if strcmpi(string(const.tracker),"NXT")==1 || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
    %pbt = table(e1name,pbins',e2name,cbins','VariableNames',["Top Pile Height","Percentage","CoTT Height","Percentage "]);
    pbt = table(e2name,cbins','VariableNames',["Pile Reveal","Percentage"]);
else
    pbt = table(e1name,pbins',e2name,cbins','VariableNames',["Top Pile Height","Percentage","CoTT Height","Percentage "]);
end
%%
for i=1:length(kpv)
    if i==1
        hagh=reshape(panels.(kpv{i}).hagh(:,:,:),[],1);
    else
        hagh_h=reshape(panels.(kpv{i}).hagh(:,:,:),[],1);
        hagh=[hagh; hagh_h];
    end
end

h1=histogram(hagh,edges,'Normalization','probability');
yl=(get(gca,'YLim'));
set(gca,'YLim',[yl(1),yl(2)+.13*yl(2)]);
ylabel('Percentage of Panels')
xlabel('Panel Height Above Ground')
n=get(h1,'Values');
xn=get(h1,"BinEdges");
formatSpec = '%.2f';
barstrings=num2str(n',formatSpec);
text(xn(2:end), n, barstrings,'horizontalalignment','center','verticalalignment','bottom','Rotation',70)

%labels=num2str([[edges(1:end-1)].' [edges(2:end)].'],'%0.2f-%0.2f');
% set labels
%set(gca,'XTickLabel',labels,'XTickLabelRotation',30)

set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);

if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' height_hist.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' height_hist.png'])
    end
    foo=1;
end
%ns slope histogram
for i=1:length(kpv)
    if i==1
        slopedist=pd_out.(kpv{i}).slp_deg;
    else
        slopedist=[slopedist;pd_out.(kpv{i}).slp_deg];
    end
end
slopedist=unique(slopedist);
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 600, 400])
edges=[floor(min(slopedist)):.25:ceil(max(slopedist))];
edges2=[floor(min(slopedist))-1:1:ceil(max(slopedist))+1];
%e2name=["3-4ft";"4-5ft";"5-6ft";"6-7ft";"7-8ft";"8-9ft";"9-10ft";"10ft+"];
[pbin pbedges]=histcounts(slopedist,edges2);
pbinp=100*(pbin/(sum(pbin(:))));
formatSpec = '%.1f';
for j=1:size(pbinp,2)
    pbins{j} = num2str(pbinp(j),formatSpec);
end
%pbt = table(e2name,pbins','VariableNames',["Pile Bin","Approximate Percentage"]);
h1=histogram(slopedist,edges,'Normalization','probability');
yl=(get(gca,'YLim'));
set(gca,'YLim',[yl(1),yl(2)+.13*yl(2)]);
ylabel('Percentage of Rows')
xlabel('Row Slope - North Positive')
n=get(h1,'Values');
xn=get(h1,"BinEdges");
formatSpec = '%.2f';
barstrings=num2str(n',formatSpec);
text(xn(2:end), n, barstrings,'horizontalalignment','center','verticalalignment','bottom','Rotation',70)
%labels=num2str([[edges(1:end-1)].' [edges(2:end)].'],'%0.2f-%0.2f');
% set labels
%set(gca,'XTickLabel',labels,'XTickLabelRotation',30)

set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);

if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' ns_hist.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' ns_hist.png'])
    end
    foo=1;
end

import mlreportgen.ppt.*
t=const.t{1};
% Create Slides
if const.simple_out==1
    if contains(string(const.tracker),"ATI",'IgnoreCase',true)
        ppt = Presentation([const.outpath '/figures/'  char(const.customer) '_' char(const.project)...
            '_Eng_RPCS_Planar_' char(t) '.pptx'],'PPT/RPCS Planar Template - ATI.potx');
    elseif matches(string(const.tracker),"NEV") || contains(string(const.project),"genpile","IgnoreCase",true)==1
        ppt = Presentation([const.outpath '/figures/'  char(const.customer) '_' char(const.project)...
            '_Eng_Study_' char(t) '.pptx'],'PPT/PPT Template - NEV.potx');
    elseif contains(string(const.tracker),"NXT","IgnoreCase",true) || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1 
        ppt = Presentation([const.outpath '/figures/'  char(const.customer) '_' char(const.project)...
            '_Eng_RPCS_Planar_' char(t) '.pptx'],'PPT/RPCS Planar Template - NXT.potx');
    end
end

open(ppt);
if const.simple_out==1
    PMN=Picture([const.outpath '/figures/' char(const.sect_name) ' height_plot.png']);
    HST=Picture([const.outpath '/figures/' char(const.sect_name) ' height_hist.png']);
    CLB=Picture([const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png']);
else
    PMN=Picture([const.outpath '/figures/' char(const.sect_name) ' height_plot.png']);
    HST=Picture([const.outpath '/figures/' char(const.sect_name) ' height_hist.png']);
    CLB=Picture([const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png']);
end
FLP=Picture('PPT/flipex.png');
if length(flipex.tpxc)>0
    replace(ppt,'Flipex',FLP)
end
FLPD=Picture('PPT/fliped.png');
if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
if ~isempty(fliped.tpxc)
    replace(ppt,'Fliped',FLPD)
end
end
replace(ppt,'PlanarMain',PMN)
replace(ppt,'Histogram',HST)
replace(ppt,'ColorBar',CLB)
%max min slope
formatSpec = ' APPROXIMATE MAX NORTH ROW SLOPE = %.1f%c \n APPROXIMATE MAX SOUTH ROW SLOPE = %.1f%c ';
northslope=max(slp_deg_ru(slp_deg_ru>0));
if isempty(northslope)==1
    northslope=0;
end
southslope=-min(slp_deg_ru(slp_deg_ru<0));
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
    %tthstrh=sprintf(formatSpec, max(pile_reveal_ru), min(pile_reveal_ru));
    tthstrh=sprintf(formatSpec, max(cott_ft_ru), min(cott_ft_ru));
else
    formatSpec = ' MAX TORQUE TUBE HEIGHT = %.1f ft \n MIN TORQUE TUBE HEIGHT = %.1f ft';
    tthstrh=sprintf(formatSpec, max(cott_ft_ru), min(cott_ft_ru));
end
tthstr = Paragraph(tthstrh);
tthstr.Font="Hack";
tthstr.FontSize="8";
replace(ppt,"MinMaxRow",tthstr)
%project name
formatSpec = 'PROJECT\n%s';
pj=upper(string(const.sect_name));
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
pause(2);
if matches(string(const.tracker),"NEV") || contains(string(const.project),"genpile","IgnoreCase",true)==1
    fname=[const.outpath '/figures/'  char(const.customer) '_' char(const.project) '_Eng_Study_' char(t) '.pptx'];
else
    fname=[const.outpath '/figures/'  char(const.customer) '_' char(const.project) '_Eng_RPCS_Planar_' char(t) '.pptx'];
    % else
    %     fname=[const.outpath '/figures/'  char(const.customer) '_' char(const.project) '_Eng_RPCS_Planar_' char(t) '.pptx'];
end

if const.simple_out==1
    delete([const.outpath '/figures/' char(const.sect_name) ' height_plot.png']);
    delete([const.outpath '/figures/' char(const.sect_name) ' height_hist.png']);
    delete([const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png']);
    % else
    %     PMN=Picture([const.outpath '/figures/' char(const.sect_name) ' height_plot.png']);
    %     HST=Picture([const.outpath '/figures/' char(const.sect_name) ' height_hist.png']);
    %     CLB=Picture([const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png']);
end
if ispc
    pptview(fname,'converttopdf')
    fnamepdf=strrep(fname,'pptx','pdf');
else
    fnamepdf=fname;
end
end



