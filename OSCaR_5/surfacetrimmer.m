function [surface]=surfacetrimmer(kpv, ...
    dat,span,const,surface)

xqc=reshape(surface.xq,[],1); %for use in delaunay
yqc=reshape(surface.yq,[],1); %for use in delaunay
%eliminate areas with no topo
indx=0;
for i=1:length(kpv)
    for j=1:length(dat.(kpv{i}).kpn)
        for k=1:length(span.(dat.(kpv{i}).span{j}).kps)
            indx=indx+1;
            xp(indx,1)=dat.(kpv{i}).tpx(j);
            yp(indx,1)=dat.(kpv{i}).tpy(j)-span.(dat.(kpv{i}).span{j}).kps(k);
        end
    end
end

DT = delaunayTriangulation(xp,yp);
[vi, d] = nearestNeighbor(DT,xqc,yqc);
surface.mask = d>50;
surface.zog=surface.zog;
surface.zog(surface.mask)=NaN;
surface.zog=surface.zog;
end