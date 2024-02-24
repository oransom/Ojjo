function plotpile_flood(const, surface, kpv, pd_out)
if strcmpi(const.flood,'na')==0
box_x=3;
box_y=170;
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 2400])
levels=min(min(surface.zog(:,:))):5:max(max(surface.zog(:,:))); %5' contours
contour3(surface.xq,surface.yq,surface.zog,levels)
foo=surface.flood;
foo(foo<0.25)=NaN;
ax1=axes;
surf(ax1, surface.xq,surface.yq,foo);
shading interp
alpha 0.3
view(2)
% hold on

for i=1:length(kpv)
    ax2=axes;
    rbx=reshape(pd_out.(kpv{i}).bpx,[],1);
    rby=reshape(pd_out.(kpv{i}).bpy,[],1);
    rbz=reshape(pd_out.(kpv{i}).bpz,[],1);
    rtx=reshape(pd_out.(kpv{i}).tpx,[],1);
    rty=reshape(pd_out.(kpv{i}).tpy,[],1);
    rtz=reshape(pd_out.(kpv{i}).tpz,[],1);
    ts=reshape(pd_out.(kpv{i}).total_steel,[],1);
    psize=length(rtx)*.1;
    scatter3(ax2,rtx,rty,ts,psize,ts,'filled');
    view(2)
end

linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'jet')
colormap(ax2,'cool')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0175 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0175 .815]);
daspect([1 1 1])

title(sprintf('Pile layout for run %s', string(const.sect_name)))
if const.writefiles == 1
    if const.simple_out==0
    mkdir([const.outpath '/figures' ]);
    saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' pile_layout_flood.png'])
    else
    mkdir([const.outpath '/figures' ]);
    saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' pile_layout_flood.png'])
    end
end
end
end
