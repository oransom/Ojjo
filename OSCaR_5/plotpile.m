function plotpile(const, surface, kpv, po, cpl, flipex, ns_ugly, pd_out)
% box_x=3;
% box_y=170;
figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 2400])
levels=min(min(surface.zog(:,:))):5:max(max(surface.zog(:,:))); %1' contours
c1=contour3(surface.xq,surface.yq,surface.zog,levels,'-k');
alpha(0.3)
hold on

for i=1:length(kpv)
    zfl=cpl.(kpv{i}).tpzf;
    zfl(isnan(zfl))=0; %this block gets rid of extra nans
    T=zfl~=0;
    n=sum(T,1);
    m=max(n);
    zfl2=nan(m,size(zfl,2));
    for k=1:size(zfl,2)
        zfl2(1:n(k),k) = zfl(T(:,k),k);
    end
    clear zfl
    zfl=zfl2';
    B=~isnan(zfl);
    Indices=arrayfun(@(x) find(B(x,:), 1, 'last'), 1:size(zfl,1));

    for j=1:size(cpl.(kpv{i}).tpx,2)
        xrp(j,1)=cpl.(kpv{i}).tpx(1,j);
        xrp(j,2)=cpl.(kpv{i}).tpx(Indices(j),j);
        yrp(j,1)=cpl.(kpv{i}).tpy(1,j);
        yrp(j,2)=cpl.(kpv{i}).tpy(Indices(j),j);
        zrp(j,1)=cpl.(kpv{i}).tpzf(1,j);
        zrp(j,2)=cpl.(kpv{i}).tpzf(Indices(j),j);
    end

    rbx=reshape(cpl.(kpv{i}).bpx,[],1);
    rby=reshape(cpl.(kpv{i}).bpy,[],1);
    rbz=reshape(cpl.(kpv{i}).bpz,[],1);
    rtx=reshape(cpl.(kpv{i}).tpx,[],1);
    rty=reshape(cpl.(kpv{i}).tpy,[],1);
    rtz=reshape(cpl.(kpv{i}).tpzf,[],1);
    mx=reshape(cpl.(kpv{i}).tbxd,[],1);
    my=reshape(cpl.(kpv{i}).tbyd,[],1);
    rbx=rmmissing(rbx);
    rby=rmmissing(rby);
    rbz=rmmissing(rbz);
    rtx=rmmissing(rtx);
    rty=rmmissing(rty);
    rtz=rmmissing(rtz);
    rtz_diff=rtz-rbz;
    mx=rmmissing(mx);
    my=rmmissing(my);
    psize=length(rtz)*.005;
    scatter3(rtx,rty,rtz,psize,[.3 .3 .3],'filled')
    colormap(jet);
    %colorbar;
    %quiver(rbx,rby,mx,my*100)
end
box_x=const.xpdim/304.8+.5;
for j=1:length(flipex.tpxc)
    box_y=flipex.row_len(j)/2+2;
    P = [flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j); flipex.tpxc(j)-box_x flipex.tpyc(j)+box_y flipex.tpzc(j);...
        flipex.tpxc(j)+box_x flipex.tpyc(j)+box_y flipex.tpzc(j); flipex.tpxc(j)+box_x flipex.tpyc(j)-box_y flipex.tpzc(j);...
        flipex.tpxc(j)-box_x flipex.tpyc(j)-box_y flipex.tpzc(j)];
    fill3(P(:,1),P(:,2),P(:,3),'r--', 'linewidth',1.5)
    %set(P,'FaceAlpha',0.3);
end
% for i=1:length(flipex.tpxc)
%     P = [flipex.tpxc(i)-box_x flipex.tpyc(i)-box_y flipex.tpzc(i); flipex.tpxc(i)-box_x flipex.tpyc(i)+box_y flipex.tpzc(i);...
%         flipex.tpxc(i)+box_x flipex.tpyc(i)+box_y flipex.tpzc(i); flipex.tpxc(i)+box_x flipex.tpyc(i)-box_y flipex.tpzc(i);...
%         flipex.tpxc(i)-box_x flipex.tpyc(i)-box_y flipex.tpzc(i)];
%     plot3(P(:,1),P(:,2),P(:,3),'b', 'linewidth',1.5)
% end


scatter3(ns_ugly.tpx,ns_ugly.tpy,ns_ugly.tpz,3*psize,'dr','filled')
view(2)
title(sprintf('Pile layout for run %s', string(const.sect_name)))

daspect([1 1 1])
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' pile_layout.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' pile_layout.png'])
    end
    foo=1;
end
