function nev_plot(kpv,pd_out,F_og,surface,ppp_out,const)
%change some variable names for compactness - not best for memory, but
%whatever
p=ppp_out; %tpx=zeros(200,200);

for i=1:length(kpv)
    edges2=[0,1.5,2.5,3.5,1000];
    e2name=["Straight";"Single";"Double";"Needs Attention"];
    [pbin pbedges]=histcounts(p.(kpv{i}).NevJoint,edges2);
    pbt = table(e2name,pbin','VariableNames',["Joint Type","Approximate Count"]);
    if const.writefiles == 1
        if const.simple_out==0
            mkdir([const.outpath '/nevados' ]);
            writetable(pbt,[const.outpath '/nevados/' char(const.sect_name) '_joint_table.csv'])
        else
            mkdir([const.outpath '/nevados' ]);
            writetable(pbt,[const.outpath '/nevados/' char(const.sect_name) '_joint_table.csv'])
        end
    end
end

 figure ('Visible', 'off')
 set(gcf, 'Position',  [100, 100, 1600, 1200])
 levels=min(min(surface.zog(:,:)))/1000000:1/1000000:max(max(surface.zog(:,:)))/1000000; %1' contours /1000000 so they are under the rows
 contour3(surface.xq,surface.yq,surface.zog/1000000,levels,':k')
 alpha(0.3)
 hold on

for i=1:length(kpv)   
        scatter3(p.(kpv{i}).tpx,p.(kpv{i}).tpy,p.(kpv{i}).NevJoint,15,p.(kpv{i}).NevJoint,'filled')
        for j=1:length(p.(kpv{i}).tpx)
            if p.(kpv{i}).NevJoint(j)==999
                scatter3(p.(kpv{i}).tpx(j),p.(kpv{i}).tpy(j),p.(kpv{i}).NevJoint(j),15,'w','filled')
            end
        end
        view(2)
        colormap jet
        %alpha(0.6)
        hold on
    minx(i)=min(p.(kpv{i}).tpx(:));
    maxx(i)=max(p.(kpv{i}).tpx(:));
    miny(i)=min(p.(kpv{i}).tpy(:));
    maxy(i)=max(p.(kpv{i}).tpy(:));
end
colorbar
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
clim([1 3])
if const.writefiles == 1
         if const.simple_out==0
             mkdir([const.outpath '/figures' ]);
             saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' nev_plot.png'])
         else
             mkdir([const.outpath '/figures' ]);
             saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' nev_plot.png'])
         end
end

%%%%%-------------------------------------XTR works? 
 figure ('Visible', 'off')
 set(gcf, 'Position',  [100, 100, 1600, 1200])
 levels=min(min(surface.zog(:,:)))/1000000:1/1000000:max(max(surface.zog(:,:)))/1000000; %1' contours /1000000 so they are under the rows
 contour3(surface.xq,surface.yq,surface.zog/1000000,levels,':k')
 alpha(0.3)
 hold on
% [ind idx]=find(pd_out.(kpv{i}).wouldXTRwork==1);
% [ind2 idx2]=find(pd_out.(kpv{i}).wouldXTRwork==0);
for i=1:length(kpv)   
        % scatter3(pd_out.(kpv{i}).bpx(ind),pd_out.(kpv{i}).bpy(ind),pd_out.(kpv{i}).wouldXTRwork(ind),5,'go','filled')
        % hold on
        % scatter3(pd_out.(kpv{i}).bpx(ind2),pd_out.(kpv{i}).bpy(ind2),pd_out.(kpv{i}).wouldXTRwork(ind2),5,'ro','filled')
        % view(2)
        % %alpha(0.6)
    minx(i)=min(p.(kpv{i}).tpx(:));
    maxx(i)=max(p.(kpv{i}).tpx(:));
    miny(i)=min(p.(kpv{i}).tpy(:));
    maxy(i)=max(p.(kpv{i}).tpy(:));
end
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
if const.writefiles == 1
         if const.simple_out==0
             mkdir([const.outpath '/figures' ]);
             saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' xtr_works.png'])
         else
             mkdir([const.outpath '/figures' ]);
             saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' xtr_works.png'])
         end
end
%%%%%-------------------------------------Nevados works? 

 figure ('Visible', 'off')
 set(gcf, 'Position',  [100, 100, 1600, 1200])
 levels=min(min(surface.zog(:,:)))/1000000:1/1000000:max(max(surface.zog(:,:)))/1000000; %1' contours /1000000 so they are under the rows
 contour3(surface.xq,surface.yq,surface.zog/1000000,levels,':k')
 alpha(0.3)
 hold on
[ind idx]=find(p.(kpv{i}).NevJoint<=3);
[ind2 idx2]=find(p.(kpv{i}).NevJoint==999);
for i=1:length(kpv)   
        scatter3(p.(kpv{i}).tpx(ind),p.(kpv{i}).tpy(ind),p.(kpv{i}).NevJoint(ind),5,'go','filled')
        hold on
        scatter3(p.(kpv{i}).tpx(ind2),p.(kpv{i}).tpy(ind2),p.(kpv{i}).NevJoint(ind2),5,'ro','filled')
        view(2)
        %alpha(0.6)
    minx(i)=min(p.(kpv{i}).tpx(:));
    maxx(i)=max(p.(kpv{i}).tpx(:));
    miny(i)=min(p.(kpv{i}).tpy(:));
    maxy(i)=max(p.(kpv{i}).tpy(:));
end
xlim([min(minx)-50 max(maxx)+50])
ylim([min(miny)-50 max(maxy)+50])
if const.writefiles == 1
         if const.simple_out==0
             mkdir([const.outpath '/figures' ]);
             saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' nevados_works.png'])
         else
             mkdir([const.outpath '/figures' ]);
             saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' nevados_works.png'])
         end
end

end