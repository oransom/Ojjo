function [sname, snamepdf]=rpcs_plot_slope_new(kpv,pd_out,F_og,surface,ppp_out,const,grd_ext,fliped, flipex)
%change some variable names for compactness - not best for memory, but
%whatever
p=ppp_out; %tpx=zeros(200,200);
dat_bins=const.cott_bins{1};
%create RGB's for dat_bins
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

yi=100;
xi=11; %should really be an odd number
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

figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 1800])
% if (max(max(surface.zog))-min(min(surface.zog)))<100
%     levels=10000+min(min(surface.zog))/1000000:1/1000000:10000+max(max(surface.zog))/1000000; %1' contours /1000000 so they are under the rows
% elseif (max(max(surface.zog))-min(min(surface.zog)))<200
%     levels=10000+min(min(surface.zog))/1000000:2/1000000:10000+max(max(surface.zog))/1000000;
% else
%     levels=10000+min(min(surface.zog))/1000000:5/1000000:10000+max(max(surface.zog))/1000000;
% end
% contour3(surface.xq,surface.yq,surface.zog/1000000,levels,':k')
% alpha(0.4)
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
            slp{:,l}=pd_out.(kpv{i}).slp_deg(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
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
            hag=eh(:,:,j)-F_og(xh(:,:,j),yh(:,:,j)); %for future FYI we can use histcounts here with edges to bin distance above ground
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
    scatter3(X,Y,Z,5,Zcol,'MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15);

    % %GRADING EXTENT PLOTTER
    % if ~isempty(grd_ext)
    %     for j=1:length(grd_ext.gbx)
    %         plot(grd_ext.xgh{j}(grd_ext.gbx{j}),grd_ext.ygh{j}(grd_ext.gby{j}),'k--','LineWidth',0.5)
    %         hold on
    %     end
    % end
    %slope boxes
    box_x=xpdim/2/304.8+2;
    for j=1:length(tpx)
        if ~isempty(tpx{j})
            if max(abs(cell2mat(slp(j))))<4.5
                Pz4=[min(cell2mat(tpx(j)))-box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))-box_x max(cell2mat(tpy(j)))+3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))+box_x max(cell2mat(tpy(j)))+3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))+box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))-box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)))];
                plot3(Pz4(:,1),Pz4(:,2),Pz4(:,3),'g-', 'linewidth',3.5)
            elseif max(abs(cell2mat(slp(j))))>4.5 && max(abs(cell2mat(slp(j))))<8.5
                P48=[min(cell2mat(tpx(j)))-box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))-box_x max(cell2mat(tpy(j)))+3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))+box_x max(cell2mat(tpy(j)))+3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))+box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))-box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)))];
                plot3(P48(:,1),P48(:,2),P48(:,3),'Color', [1 .6412 0], 'linewidth',3.5)
            else
                P8P=[min(cell2mat(tpx(j)))-box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))-box_x max(cell2mat(tpy(j)))+3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))+box_x max(cell2mat(tpy(j)))+3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))+box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)));...
                    min(cell2mat(tpx(j)))-box_x min(cell2mat(tpy(j)))-3 max(cell2mat(top(j)))];
                plot3(P8P(:,1),P8P(:,2),P8P(:,3),'r-', 'linewidth',3.5)
            end

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
    if const.min_wp ~= const.max_wp
        clim([const.min_wp const.max_wp])
    else
        clim([const.min_wp+offset-.5 const.max_wp+offset+.5])
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
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
if const.min_wp ~= const.max_wp
    clim([const.min_wp const.max_wp])
else
    clim([const.min_wp+offset-.5 const.max_wp+offset+.5])
end

% %flipex boxes
% box_x=xpdim/304.8+.5;
% if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
%     for j=1:length(fliped.tpxc)
%         box_y=fliped.row_len(j)/2+2;
%         P = [fliped.tpxc(j)-box_x fliped.tpyc(j)-box_y fliped.tpzc(j); fliped.tpxc(j)-box_x fliped.tpyc(j)+box_y fliped.tpzc(j);...
%             fliped.tpxc(j)+box_x fliped.tpyc(j)+box_y fliped.tpzc(j); fliped.tpxc(j)+box_x fliped.tpyc(j)-box_y fliped.tpzc(j);...
%             fliped.tpxc(j)-box_x fliped.tpyc(j)-box_y fliped.tpzc(j)];
%         plot3(P(:,1),P(:,2),P(:,3),'r--', 'linewidth',1.5)
%     end
%     for j=1:length(flipex.tpxc)
%         box_y=flipex.row_len(j)/2+2;
%         P = [flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j); flipex.tpxc(j)-box_x flipex.tpyc(j)+box_y flipex.tpzc(j);...
%             flipex.tpxc(j)+box_x flipex.tpyc(j)+box_y flipex.tpzc(j); flipex.tpxc(j)+box_x flipex.tpyc(j)-box_y flipex.tpzc(j);...
%             flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j)];
%         plot3(P(:,1),P(:,2),P(:,3),'r-', 'linewidth',1.5)
%     end
% else
%     for j=1:length(flipex.tpxc)
%         box_y=flipex.row_len(j)/2+2;
%         P = [flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j); flipex.tpxc(j)-box_x flipex.tpyc(j)+box_y flipex.tpzc(j);...
%             flipex.tpxc(j)+box_x flipex.tpyc(j)+box_y flipex.tpzc(j); flipex.tpxc(j)+box_x flipex.tpyc(j)-box_y flipex.tpzc(j);...
%             flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j)];
%         plot3(P(:,1),P(:,2),P(:,3),'r-', 'linewidth',1.5)
%     end
% end

if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' rowslope_plot.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' rowslope_plot.png'])
    end
end

%% Other Plots
% %Colorbar for height above ground plot
% fig_leg=figure('Visible', 'off'); %create stand alone color bar
% axis off
% for k=1:length(dat_bins.bin_start)
%     fig_leg=scatter(nan,nan, 'MarkerFaceColor', dat_bins.rgb{k},'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
%     hold on
% end
% legend(dat_bins.bin_name')
% % Initial values to capture the entire legend
% % Should fit most modern screens
% set(gcf,'Position',[0,0,1024,1024]);
% % Call the legend to your choice, I used a horizontal legend here
% legend_handle = legend('Orientation','vertical');
% % Set the figure Position using the normalized legend Position vector
% % as a multiplier to the figure's current position in pixels
% % This sets the figure to have the same size as the legend
% set(gcf,'Position',(get(legend_handle,'Position')...
%     .*[0, 0, 1, 1].*get(gcf,'Position')));
% % The legend is still offset so set its normalized position vector to
% % fill the figure
% set(legend_handle,'Position',[0,0,1,1]);
% % Put the figure back in the middle screen area
% set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);
% 
% if const.writefiles == 1
%     if const.simple_out==0
%         mkdir([char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/general/' char(const.sect_name) '/figures' ]);
%         saveas(gcf,[char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/general/' char(const.sect_name) '/figures/' char(const.sect_name) ' height_colorbar.png'])
%     else
%         mkdir([char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/' char(const.sect_name) '/figures' ]);
%         saveas(gcf,[char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/' char(const.sect_name) '/figures/' char(const.sect_name) ' height_colorbar.png'])
%     end
% end

formatSpec = '%.1f';
for j=1:length(slp)
    if length(slp{:,j})>1
        slpm(j)=max(abs(cell2mat(slp(:,j))));
    else
        slpm(j)=NaN;
    end
end
sl4=100*(sum(slpm<4.5)/length(slpm));
s4t8=100*(sum(slpm>4.5 & slpm<8.5)/length(slpm));
sg8=100*(sum(slpm>8.5)/length(slpm));
tcut = num2str(sl4,formatSpec);
tfill = num2str(s4t8,formatSpec);
tarea = num2str(sg8,formatSpec);

import mlreportgen.ppt.*
t=const.t{1};
pbt = table(string(tcut),string(tfill),string(tarea),'VariableNames',["Rows 0˚-4.5˚","Rows 4.5˚-8.5˚","Rows > 8.5˚"]);
ppt = Presentation([const.outpath '/figures/'  char(const.customer) '_' char(const.project)...
            '_ROW_SLOPE_' char(t) '.pptx'],'PPT/rowslope_fig.potx');
open(ppt);
PMN=Picture([const.outpath '/figures/' char(const.sect_name) ' rowslope_plot.png']);
%CLB=Picture([char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/' char(const.sect_name) '/figures/' char(const.sect_name) ' height_colorbar.png']);
%nsHST=Picture([char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/' char(const.sect_name) '/figures/' char(const.sect_name) ' ns_hist.png']);
%ewHST=Picture([char(const.outpath{1}) '/' char(const.project{1}) '_output_' char(const.t{1}) '/' char(const.sect_name) '/figures/' char(const.sect_name) ' ew_hist.png']);
replace(ppt,'PlanarMain',PMN)
%replace(ppt,'ColorBar',CLB)
%replace(ppt,'EWHistogram',ewHST)
%replace(ppt,'NSHistogram',nsHST)
% FLP=Picture('PPT/flipex.png');
% if length(flipex.tpxc)>0
%     replace(ppt,'Flipex',FLP)
% end
% FLPD=Picture('PPT/fliped.png');
% if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
% if ~isempty(fliped.tpxc)
%     replace(ppt,'Fliped',FLPD)
% end
% end
%max min slope
formatSpec = ' APPROXIMATE MAX NORTH ROW SLOPE = %.1f%c \n APPROXIMATE MAX SOUTH ROW SLOPE = %.1f%c ';
slpstrh=sprintf(formatSpec, max(slp_deg_ru), char(176), -min(slp_deg_ru), char(176));
slpstr = Paragraph(slpstrh);
slpstr.Font="Hack";
slpstr.FontSize="8";
replace(ppt,"RowSlope",slpstr)
%max min tth
formatSpec = ' MAX TORQUE TUBE HEIGHT = %.1f ft \n MIN TORQUE TUBE HEIGHT = %.1f ft';
tthstrh=sprintf(formatSpec, max(cott_ft_ru), min(cott_ft_ru));
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


sname=([const.outpath '/figures/'  char(const.customer) '_' char(const.project) '_ROW_SLOPE_' char(t) '.pptx']);

if ispc
    pptview(sname,'converttopdf')
    snamepdf=strrep(sname,'pptx','pdf');
else
    snamepdf=sname;
end
delete([const.outpath '/figures/' char(const.sect_name) ' height_colorbar.png']);
delete([const.outpath '/figures/' char(const.sect_name) ' rowslope_plot.png']);
delete([const.outpath '/figures/' char(const.sect_name) ' ns_hist.png']);
delete([const.outpath '/figures/' char(const.sect_name) ' ew_hist.png']);