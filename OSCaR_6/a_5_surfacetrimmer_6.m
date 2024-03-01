function [surface]=a_5_surfacetrimmer_6(dpile,surface)
xqc=reshape(surface.xq,[],1); %for use in delaunay
yqc=reshape(surface.yq,[],1); %for use in delaunay
DT = delaunayTriangulation(dpile.bpx,dpile.bpy);
[~, d] = nearestNeighbor(DT,xqc,yqc);
surface.mask = d>100;
surface.zog(surface.mask)=NaN;
end