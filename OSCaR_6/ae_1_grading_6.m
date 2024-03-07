
function [surface, grading] = ae_1_grading_6(const,surface,dpileh)

surface.graded=NaN(size(surface.zog));
padsizex=0.3*const.xspacing;
padsizey=15;

cut=find(dpileh.cotth<const.min_wp); %table index for cut
fill=find(dpileh.cotth>const.max_wp); %table index for fill

cuth=dpileh.bpzc(cut)+(dpileh.cotth(cut)-const.min_wp);
fillh=dpileh.bpzc(fill)+(dpileh.cotth(fill)-const.max_wp);

idcf=[cut;fill];
amcf=[cuth;fillh];

if isempty(idcf) %no grading required
    fprintf('no grading \n')
    grading.cut=0;
    grading.fill=0;
    grading.ext=[];
    grading.area=0;
    grading.xyz=[];
    grading.gname=[];
    grading.gnamepdf=[];
    surface.graded=surface.zog;
    surface.Fg=surface.F_og;
    surface.sdiff=surface.graded-surface.zog;
    surface.sdiff(abs(surface.sdiff)<0.01)=NaN;
    return
end
%grading areas
gext_xm=dpileh.tpxc(idcf)-padsizex;
gext_xp=dpileh.tpxc(idcf)+padsizex;
gext_ym=dpileh.tpyc(idcf)-padsizey;
gext_yp=dpileh.tpyc(idcf)+padsizey;
%padded grading areas
gext_xmp=dpileh.tpxc(idcf)-2*padsizex;
gext_xpp=dpileh.tpxc(idcf)+2*padsizex;
gext_ymp=dpileh.tpyc(idcf)-2*padsizey;
gext_ypp=dpileh.tpyc(idcf)+2*padsizey;
for i=1:height(idcf)
    pgrad_areas=find(surface.xq>gext_xmp(i) & surface.xq<gext_xpp(i) & surface.yq>gext_ymp(i) & surface.yq<gext_ypp(i)); %padded area
    surface.graded(pgrad_areas)=surface.zog(pgrad_areas);
end
for i=1:height(idcf)
    grad_areas=find(surface.xq>gext_xm(i) & surface.xq<gext_xp(i) & surface.yq>gext_ym(i) & surface.yq<gext_yp(i)); %strict grading area
    surface.graded(grad_areas)=amcf(i);
end
surface.graded=imfilter(surface.graded,fspecial('average',[const.gsx const.gsy]),"replicate");
surface.graded(isnan(surface.graded))=surface.zog(isnan(surface.graded));
surface.sdiff=surface.graded-surface.zog;
surface.sdiff(abs(surface.sdiff)<0.01)=NaN;
x=reshape(surface.xq,[],1);
y=reshape(surface.yq,[],1);
z=reshape(surface.zog,[],1);
zg=reshape(surface.graded,[],1);
zd=reshape(surface.sdiff,[],1);

xyzo=[x,y,z]; xyzo=rmmissing(xyzo);
xyzd=[x,y,zd]; xyzd=rmmissing(xyzd);
xyzg=[x,y,zg]; xyzg=rmmissing(xyzg);

if const.writesurf ==1
    writematrix(xyzd, append(const.gpath{1}, '/', 'graded_areas.csv'),'delimiter',',')
    writematrix(xyzg, append(const.gpath{1}, '/', 'fullsurf.csv'),'delimiter',',')
    writematrix(xyzo, append(const.gpath{1}, '/', 'eg_mat_surf.csv'),'delimiter',',')
end

grading.fill=sum(sum(surface.sdiff(surface.sdiff > 0)))*(const.bin_x*const.bin_y)/27; %fill
grading.cut=abs(sum(sum(surface.sdiff(surface.sdiff < 0)))*(const.bin_x*const.bin_y)/27); %cut
grading.area=(sum(sum(surface.sdiff>0))+abs(sum(sum(surface.sdiff<0))))*(const.bin_x*const.bin_y)/43560;

surface.Fg=scatteredInterpolant(x,y,zg,'natural','none');
dpileh.cottpg_check=dpileh.tpzc-surface.Fg(dpileh.tpxc,dpileh.tpyc);

end