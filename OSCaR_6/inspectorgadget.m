
cott=dpile.cott_pg;
flddep=surface.Flood_s(dpile.tpxc,dpile.tpyc)-surface.F_og(dpile.tpxc,dpile.tpyc);
cott2fld=cott-flddep;
figure visible on
scatter(find(cott),sort(cott));
title('POST GRADE COTT HEIGHTS');
% figure
% scatter(find(cott2fld),sort(cott2fld));
fldid=find(cott2fld<const.freeboard);
figure
scatter3(dpile.tpxc(fldid),dpile.tpyc(fldid),dpile.tpzc(fldid)+1000,5,"black","filled");
hold on
redChannel = surface.sdiff < 0; %this is cut
greenChannel = surface.sdiff > 0; %this is fill
blueChannel = zeros(size(greenChannel));
colors = double(cat(3, redChannel, greenChannel, blueChannel));
s=surf(surface.xq, surface.yq, surface.sdiff, colors);
s.AlphaData = (surface.sdiff<=-0.001)|(surface.sdiff>=0.001);
s.FaceAlpha = 0.35;
s.EdgeAlpha = 0.5;
daspect([1 1 1])
shading interp
view(2)

badfld=flddep(fldid);
grdlvf=surface.Fg(dpile.tpxc(fldid),dpile.tpyc(fldid));
fldlv=surface.Flood_s(dpile.tpxc(fldid),dpile.tpyc(fldid));

cott_flood=dpile.tpzc-surface.Flood_s(dpile.tpxc,dpile.tpyc);
figure
scatter(find(cott_flood),sort(cott_flood));
title('COTT ABOVE FLOOD FINAL PG HEIGHTS')
for i=1:height(drow)
    if ~isnan(drow.nne(i))
        re(i)=abs(drow.rowzavg(i)-drow.rowzavg(drow.nne(i)));
    else
        re(i)=NaN;
    end
    if ~isnan(drow.nnw(i))
        rw(i)=abs(drow.rowzavg(i)-drow.rowzavg(drow.nnw(i)));
    else
        rw(i)=NaN;
    end
end

figure
scatter(find(re),sort(re));
title('ROW EAST DELTA')
figure
scatter(find(rw),sort(rw));
title('ROW WEST DELTA')