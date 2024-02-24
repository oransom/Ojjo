function [dat] = NPtrimmer(kpv,dat_in,span,const,surface)

for i=1:length(kpv)
    datin=dat_in.(kpv{i});
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
    for j=1:length(datin.span)
        datin.row_length(j)=span.(datin.span{j}).kps(end);
    end
    datin.in = inpolygon(datin.tpx,datin.tpy,foox(k),fooy(k)) ...
        & inpolygon(datin.tpx,datin.tpy-datin.row_length,foox(k),fooy(k));
    intdat=double(datin.in);
    out=nnz(~intdat);
    fprintf('%5.0f NP fall outside topo boundary for project %s section %s \n',out,char(const.sect_name),(kpv{i}));
    intdat(intdat==0)=NaN;
    datin.in=intdat;
    datin=rmmissing(datin);
    datin = removevars(datin,{'in'}); %remove the in column, because why not
    dat.(kpv{i})=datin;
end