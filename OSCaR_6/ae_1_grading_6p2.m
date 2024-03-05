function [surface, dpileh, grading] = ae_1_grading_6p2(const,surface,dpileh,drowh)

sgnan=NaN(size(surface.zog));
padsizex=0.3*const.xspacing;
padsizey=15;
tcfr=const.gcutfill;
%-----
cut=find(dpileh.cotth<const.min_wp); %pile table index for cut
fill=find(dpileh.cotth>const.max_wp); %pile table index for fill
g.idcf=[cut;fill];
%-----
cuta=dpileh.cotth(cut)-const.min_wp; %total cut required in feet
filla=dpileh.cotth(fill)-const.max_wp; %total fill required in feet
g.mcf=[cuta;filla];
%-----
cuth=dpileh.bpzc(cut)+(dpileh.cotth(cut)-const.min_wp);
fillh=dpileh.bpzc(fill)+(dpileh.cotth(fill)-const.max_wp);
g.amcf=[cuth;fillh];
%-----
g=struct2table(g);
g=sortrows(g);

dpileh.cottcg=dpileh.cotth;

if isempty(g.idcf) %no grading required
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

row.rrg=zeros(length(g.idcf),1); row.mag=zeros(length(g.idcf),1); %row lookup for drow, magnitude of grading
row.d2w=zeros(length(g.idcf),1); row.d2e=zeros(length(g.idcf),1); %delta to east delta to west
row.lmb=zeros(length(g.idcf),1); row.umb=zeros(length(g.idcf),1); %lower movement bound upper movement bound

for i=1:height(g)
    row.rrg(i,1)=find(g.idcf(i)>=drowh.si & g.idcf(i)<=drowh.ei); %rows requring grading
    row.mag=-g.mcf; %dumb way to do this, but making it readable - how much the row needs to move to get parity. Flipping sign so that
    %mag and allow contridict for ease of calcuation below
    if const.flipex==1
        row.flp(i,1)=drowh.flip2ext(row.rrg(i)); %how much can the row move
    else
        row.flp(i,1)=2;
    end
    %sign convention looking both west and east is important to be the same
    %so we always know if target row is above or below both rows to each
    %side + means above that direction - means below that direction
    if ~isnan(drowh.nnw(row.rrg(i))) && (drowh.ntpxc(row.rrg(i))-drowh.ntpxc(drowh.nnw(row.rrg(i))))<const.xspacing+10
        row.d2w(i,1)=drowh.rowzavg(row.rrg(i))-drowh.rowzavg(drowh.nnw(row.rrg(i)));
    end
    if ~isnan(drowh.nne(row.rrg(i))) && (drowh.ntpxc(drowh.nne(row.rrg(i)))-drowh.ntpxc(row.rrg(i)))<const.xspacing+10
        row.d2e(i,1)=drowh.rowzavg(row.rrg(i))-drowh.rowzavg(drowh.nne(row.rrg(i)));
    end
    hghr=max(row.d2w(i),row.d2e(i));
    lwr=min(row.d2w(i),row.d2e(i));
    row.lmb(i)=hghr-const.row_delta;
    row.umb(i)=lwr+const.row_delta;
end
row=struct2table(row);
idcfu=unique(row.rrg); %unique rows from row table that require grading
for i=1:numel(idcfu)
    j=idcfu(i); %shorthand
    tar.row(i,1)=(j);
    tar.min(i,1)=const.min_wp;
    tar.max(i,1)=const.max_wp;
    tar.ravg(i,1)=mean(dpileh.cotth(drowh.si(j):drowh.ei(j)));
    tar.rvavg(i,1)=mean(row.mag(row.rrg==idcfu(i)));
    tar.lmb(i,1)=mean(row.lmb(row.rrg==idcfu(i)));
    tar.umb(i,1)=mean(row.umb(row.rrg==idcfu(i)));
    if ~strcmpi(const.flood,'na')
        tar.fmtth(i,1)=min(dpileh.cotth(j)-dpileh.matpzc(j));
    else
        tar.fmtth(i,1)=tar.lmb(i,1);
    end
    tar.lmb(i)=min(tar.lmb(i),tar.fmtth(i));
    %if rvavg is positive move it down to increase cut and lower fill
    %if rvavg is negative move it up to increase fill and lower cut
    %rvavg of zero is balanced cut fill
    if tar.rvavg>0
        tar.tmove(i,1)=max(-tar.rvavg(i),tar.lmb(i));
    else
        tar.tmove(i,1)=min(-tar.rvavg(i),tar.umb(i));
    end
end
tar=struct2table(tar);
cottcgh=dpileh.cottcg;
c2fr=0.9*tcfr;%grading.cut/grading.fill - set to run at least once
cfsso=0.1; %initial cut fill step size
cfss=cfsso; %for first step
c2fdo=0; %for first step
iter=0;

while c2fr<tcfr-0.01 || c2fr>tcfr+0.01
    iter=iter+1;
    fprintf('loop %2.0f c2fr %3.2f %3.2f s \n',iter,c2fr,toc);
    surface.graded=sgnan;
    c2fro=c2fr;
    for i=1:height(tar)
        j=tar.row(i); %shorthand for the row we're operating on
        tar.tmoved(i)=tar.tmove(i);
        dpileh.cottcg(drowh.si(j):drowh.ei(j))=cottcgh(drowh.si(j):drowh.ei(j))-tar.tmove(i); %this zeros out the pile - in theory
        %also - we're always referencing the balanced cott heights
        %calculated above
    end
    %-----
    cutg=find(dpileh.cottcg<const.min_wp); %table index for cut
    fillg=find(dpileh.cottcg>const.max_wp); %table index for fill
    wg.idcf=[cutg;fillg];
    %-----
    cuta=dpileh.cottcg(cutg)-const.min_wp; %total cut required in feet
    filla=dpileh.cottcg(fillg)-const.max_wp; %total fill required in feet
    wg.mcf=[cuta;filla];
    %-----
    cuth=dpileh.bpzc(cutg)+cuta;
    fillh=dpileh.bpzc(fillg)+filla;
    wg.amcf=[cuth;fillh];
    %-----
    wg=struct2table(wg);
    wg=sortrows(wg);
    %grading areas
    gext_xm=dpileh.tpxc(wg.idcf)-padsizex;
    gext_xp=dpileh.tpxc(wg.idcf)+padsizex;
    gext_ym=dpileh.tpyc(wg.idcf)-padsizey;
    gext_yp=dpileh.tpyc(wg.idcf)+padsizey;
    %padded grading areas
    gext_xmp=dpileh.tpxc(wg.idcf)-2*padsizex;
    gext_xpp=dpileh.tpxc(wg.idcf)+2*padsizex;
    gext_ymp=dpileh.tpyc(wg.idcf)-2*padsizey;
    gext_ypp=dpileh.tpyc(wg.idcf)+2*padsizey;
    for i=1:height(wg)
        pgrad_areas=find(surface.xq>gext_xmp(i) & surface.xq<gext_xpp(i) & surface.yq>gext_ymp(i) & surface.yq<gext_ypp(i)); %padded area
        surface.graded(pgrad_areas)=surface.zog(pgrad_areas);
    end
    for i=1:height(wg)
        grad_areas=find(surface.xq>gext_xm(i) & surface.xq<gext_xp(i) & surface.yq>gext_ym(i) & surface.yq<gext_yp(i)); %strict grading area
        surface.graded(grad_areas)=wg.amcf(i);
    end
    sghold=surface.graded;
    surface.graded(isnan(surface.graded))=surface.zog(isnan(surface.graded));
    surface.sdiff=surface.graded-surface.zog;
    surface.sdiff(abs(surface.sdiff)<0.01)=NaN;

    grading.fill=sum(sum(surface.sdiff(surface.sdiff > 0)))*(const.bin_x*const.bin_y)/27; %fill
    grading.cut=abs(sum(sum(surface.sdiff(surface.sdiff < 0)))*(const.bin_x*const.bin_y)/27); %cut
    grading.area=(sum(sum(surface.sdiff>0))+abs(sum(sum(surface.sdiff<0))))*(const.bin_x*const.bin_y)/43560;
    
    
    %cut 2 fill ration (c2fr) outside allowable target cut fill ratio tcfr
    c2fr=grading.cut/grading.fill;
    c2fd=c2fr-tcfr; %cut to fill delta
    cfss=abs(c2fd*cfsso);
    if abs(c2fd-c2fdo)<0.001
        break
    end
    c2fdo=c2fd;
    for i=1:height(tar)
       if c2fd>0 %if more cut than desired
           if tar.tmove(i)-cfss>tar.lmb(i) %&& tar.tmove(i)-cfss>=0
               tar.tmove(i)=tar.tmove(i)-cfss;
           end
       else %if more fill than desired
           if tar.tmove(i)+cfss<tar.umb(i) %&& tar.tmove(i)+cfss<=0
               tar.tmove(i)=tar.tmove(i)+cfss;
           end
       end
    end
    clear gext* *grad_areas wg
end
clear surface.graded surface.sdiff
fprintf('Final Ratio: %3.2f %3.2f s, cut %6.2f, fill %6.2f  \n',c2fr,toc,grading.cut,grading.fill);

sghold=imfilter(sghold,fspecial('average',[const.gsx const.gsy]),"replicate");
sghold(isnan(sghold))=surface.zog(isnan(sghold));
surface.graded=sghold;
surface.sdiff=surface.graded-surface.zog;
surface.sdiff(abs(surface.sdiff)<0.01)=NaN;
grading.fill=sum(sum(surface.sdiff(surface.sdiff > 0)))*(const.bin_x*const.bin_y)/27; %fill
grading.cut=abs(sum(sum(surface.sdiff(surface.sdiff < 0)))*(const.bin_x*const.bin_y)/27); %cut
grading.area=(sum(sum(surface.sdiff>0))+abs(sum(sum(surface.sdiff<0))))*(const.bin_x*const.bin_y)/43560;
x=reshape(surface.xq,[],1);
y=reshape(surface.yq,[],1);
z=reshape(surface.zog,[],1);
zg=reshape(surface.graded,[],1);
zd=reshape(surface.sdiff,[],1);

xyzo=[x,y,z]; xyzo=rmmissing(xyzo);
xyzd=[x,y,zd]; xyzd=rmmissing(xyzd);
xyzg=[x,y,zg]; xyzg=rmmissing(xyzg);
surface.Fg=scatteredInterpolant(x,y,zg,'natural','none');

if const.writesurf ==1
    writematrix(xyzd, append(const.gpath{1}, '/', 'graded_areas.csv'),'delimiter',',')
    writematrix(xyzg, append(const.gpath{1}, '/', 'fullsurf.csv'),'delimiter',',')
    writematrix(xyzo, append(const.gpath{1}, '/', 'eg_mat_surf.csv'),'delimiter',',')
end
end
%-----
% for i=1:numel(sects)
%     sect.cut(i,1)=sum(cuta(dpileh.sect(cut)==sects(i))); %section cut
%     sect.fill(i,1)=sum(filla(dpileh.sect(fill)==sects(i))); %section fill
%     sect.cfr(i,1)=sect.cut(i)/sect.fill(i); %section cut/fill ratio
%     sect.ach(i,1)=mean(dpileh.cotth(dpileh.sect==sects(i))); %section average cott
%     sect.ctp(i,1)=sum(dpileh.sect(cut)==sects(i)); %number of pile causing cut
%     sect.flp(i,1)=sum(dpileh.sect(fill)==sects(i)); %number of pile causing fill
%     sect.cpp(i,1)=sect.cut(i)/sect.ctp(i); %how much cut per pile
%     sect.fpp(i,1)=sect.fill(i)/sect.flp(i); %how much fill per pile
% end
% sect=struct2table(sect);
%tr=1;
% ystep=drowh.ntpyc(tr):-const.bin_y:drowh.stpyc(tr);
% xstep=ones(1,numel(ystep)).*drowh.ntpxc(tr);
% nbpz=surface.F_og(drowh.ntpxc(tr),drowh.ntpyc(tr));
% sbpz=surface.F_og(drowh.stpxc(tr),drowh.stpyc(tr));
% zh=(nbpz-sbpz)/(numel(ystep)-1);
% zlin=nbpz:-zh:sbpz;
% zsur=surface.F_og(xstep,ystep);
% figure Visible on
% scatter(ystep,zlin);
% hold on
% scatter(ystep,zsur);
% scatter(dpileh.tpyc(drowh.si(tr):drowh.ei(tr)),dpileh.tpzc(drowh.si(tr):drowh.ei(tr)))