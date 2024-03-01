function [dpile] = a_6_pointstrimmer_6(dpile,surface);

foo.x=reshape(surface.xq,[],1);
foo.y=reshape(surface.yq,[],1);
foo.z=reshape(surface.zog,[],1);
foo_tab=struct2table(foo);
foo_real=rmmissing(foo_tab);
foox=foo_real.x;
fooy=foo_real.y;
k=boundary(foox,fooy,1);
%next line returns in indicating if the query points specified by first two
%are inside or on the edge of the polygon area defined by last two.
pbound = inpolygon(dpile.tpx,dpile.tpy,foox(k),fooy(k));
intdat=double(pbound);
out=nnz(~intdat);
fprintf('%5.0f pile points fall outside topo boundary \n',out);
intdat(intdat==0)=NaN;
dpile.inbounds=intdat;
dpile=rmmissing(dpile);