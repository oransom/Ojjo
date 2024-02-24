function [total_grad,grd_ext,grd_area,xyz,gname,gnamepdf] = bfgrade3(const, kpv, po, surface, ppp_out, F_og)
p=ppp_out;
padsizex=13;
padsizey=15;
surf_holder=surface.zog;
surf_holder_trimmed=surface.zog;
fullsurf=surf_holder_trimmed;
surfmask=zeros(size(surface.xq));

ss=0.2*const.bin_x;  %rise/run 3:1=.33, 4:1=.25, 5:1=.2
%surfzero=surfmask;


for i=1:length(kpv)
    %pile height
    pile_h=ppp_out.(kpv{i}).reveal_ht_ft;
    %determine cut
    cutx=ppp_out.(kpv{i}).tpx(pile_h<const.min_wp); 
    cuty=ppp_out.(kpv{i}).tpy(pile_h<const.min_wp);
    cutz=F_og(cutx,cuty)+(pile_h(pile_h<const.min_wp))-const.min_wp-0.1; %padding by a tenth of a foot
    %determine fill
    fillx=ppp_out.(kpv{i}).tpx(pile_h>const.max_wp);
    filly=ppp_out.(kpv{i}).tpy(pile_h>const.max_wp);
    fillz=F_og(fillx,filly)+pile_h(pile_h>const.max_wp)-const.max_wp+0.1; %padding by a tenth of a foot

    %create xyz cut and fill arrays
    if i==1
        cxa=cutx; cya=cuty; cza=cutz;
        fxa=fillx; fya=filly; fza=fillz;
    else
        cxa=[cxa;cutx]; cya=[cya;cuty]; cza=[cza;cutz];
        fxa=[fxa;fillx]; fya=[fya;filly]; fza=[fza;fillz];
    end
    clear pile_h cutx cuty cutz fillx filly fillz
end

cfx=[cxa; fxa];
cfy=[cya; fya];
cfz=[cza; fza];

if isempty(cfz)==1 %no grading required
    fprintf('no grading \n')
    total_grad.cut=0;
    total_grad.fill=0;
    grd_ext=[];
    grd_area=0;
    xyz=[];
    gname=[];
    gnamepdf=[];
else
    if isempty(cxa)==1
        surfmask1=zeros(size(surfmask));
    else
        cut_size=size(cxa,1);
        for j=1:cut_size
            if j==1
                surfmask1=cza(j)*((surface.xq>=(cxa(j)-padsizex) & surface.xq<=(cxa(j)+padsizex)...
                    & surface.yq>=(cya(j)-padsizey) & surface.yq<=(cya(j)+padsizey)));
            else
                surfmaskh=cza(j)*((surface.xq>=(cxa(j)-padsizex) & surface.xq<=(cxa(j)+padsizex)...
                    & surface.yq>=(cya(j)-padsizey) & surface.yq<=(cya(j)+padsizey)));
                surfmask1(surfmaskh==cza(j))=cza(j);
            end
        end
    end

    if isempty(fxa)==1
        surfmask2=zeros(size(surfmask));
    else
        fill_size=size(fxa,1);
        for j=1:fill_size
            if j==1
                surfmask2=fza(j)*((surface.xq>=(fxa(j)-padsizex) & surface.xq<=(fxa(j)+padsizex)...
                    & surface.yq>=(fya(j)-padsizey) & surface.yq<=(fya(j)+padsizey)));
            else
                surfmaskh=fza(j)*((surface.xq>=(fxa(j)-padsizex) & surface.xq<=(fxa(j)+padsizex)...
                    & surface.yq>=(fya(j)-padsizey) & surface.yq<=(fya(j)+padsizey)));
                surfmask2(surfmaskh==fza(j))=fza(j);
            end
        end
    end

    cut_sm1=surfmask1~=0;
    fill_sm2=surfmask2~=0;
    %this order prioritizes fill over cut
    surfmask(cut_sm1)=surfmask1(cut_sm1);
    surfmask(fill_sm2)=surfmask2(fill_sm2);
    surfmask_insert=surfmask~=0;
    surf_holder(surfmask_insert)=surfmask(surfmask_insert);%surf_holder+surfmask;


    surfnonzero=surfmask~=0;
    gradx=surface.xq(surfnonzero);
    grady=surface.yq(surfnonzero);
    DT = delaunayTriangulation(gradx,grady);
    xqc=reshape(surface.xq,[],1); %for use in delaunay
    yqc=reshape(surface.yq,[],1); %for use in delaunay
    [vi, d] = nearestNeighbor(DT,xqc,yqc);
    grade_mask = d>10;
    surf_holder_trimmed=surf_holder;
    surf_holder_trimmed(grade_mask)=NaN;
    surf_holder_trimmed=imfilter(surf_holder_trimmed,fspecial('average',[const.gsx const.gsy]),"replicate"); %somehow have to control the stamp size
    clear surfnonzero
    surfnonzero=surf_holder_trimmed>0;
    fullsurf(surfnonzero)=surf_holder_trimmed(surfnonzero);

    xo=reshape(surface.xq,[],1); yo=reshape(surface.yq,[],1); zo=reshape(surface.zog,[],1);
    x=reshape(surface.xq,[],1); y=reshape(surface.yq,[],1); z=reshape(fullsurf,[],1);%z=reshape(surf_holder_trimmed,[],1);

    z_diff=z-zo;
    z_diff(z_diff==0)=NaN;

    xyzo=[xo, yo, zo];
    xyz=[x, y, z];
    xyz_diff=[x, y, z_diff];
    xyzo=rmmissing(xyzo);
    xyz=rmmissing(xyz);
    xyz_diff=rmmissing(xyz_diff);
    %find grading extent boundaries
    [grd_ext] = grading_extents(xyz_diff);
    %grd_ext=[];

    if const.writefiles == 1
        if const.writesurf ==1
            mkdir([const.outpath '/grading/']);
            dlmwrite([const.outpath '/grading/' 'graded_areas.csv'],xyz_diff,'delimiter',',','precision','%.4f')
            dlmwrite([const.outpath '/grading/' 'fullsurf.csv'],xyz,'delimiter',',','precision','%.4f')
            dlmwrite([const.outpath '/grading/' 'eg_mat_surf.csv'],xyzo,'delimiter',',','precision','%.4f')
        end
    end

    numsect=1:max(p.(kpv{i}).section);
    for j=1:max(numsect)
        if ismember(j,p.(kpv{i}).section)
            numrow(j)=max(p.(kpv{i}).row_number(p.(kpv{i}).section==j));
        else
            numrow(j)=NaN;
        end
    end
    % l=0;
    % for j=1:max(numsect)
    %     for k=1:max(numrow(j))
    %         l=l+1;
    %         tpx{:,l}=p.(kpv{i}).tpx(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
    %         tpy{:,l}=p.(kpv{i}).tpy(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
    %         top{:,l}=p.(kpv{i}).reveal_ht_ft(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
    %         elev{:,l}=p.(kpv{i}).top_pile_elev_ft(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
    %     end
    % end
    % minx(i)=min(p.(kpv{i}).tpx(:));
    % maxx(i)=max(p.(kpv{i}).tpx(:));
    % miny(i)=min(p.(kpv{i}).tpy(:));
    % maxy(i)=max(p.(kpv{i}).tpy(:));

    minx=min(min(surface.xq));
    maxx=max(max(surface.xq));
    miny=min(min(surface.yq));
    maxy=max(max(surface.yq));

    foo=surf_holder_trimmed-surface.zog;
    foop=sum(sum(foo(foo > 0)))*(const.bin_x*const.bin_y)/27; %fill
    foom=abs(sum(sum(foo(foo < 0)))*(const.bin_x*const.bin_y)/27); %cut
    fooa=(sum(sum(foo>0))+abs(sum(sum(foo<0))))*(const.bin_x*const.bin_y)/43560;
    grd_area=fooa;

    figure ('Visible', 'off')
    set(gcf, 'Position',  [100, 100, 2400, 2400])
    redChannel = foo < 0; %this is cut
    greenChannel = foo > 0; %this is fill
    blueChannel = zeros(size(greenChannel));
    colors = double(cat(3, redChannel, greenChannel, blueChannel));
    if (max(max(fullsurf(:,:)))-min(min(fullsurf(:,:))))<100
        levels=10000+min(min(fullsurf(:,:)))/1000000:1/1000000:10000+max(max(fullsurf(:,:)))/1000000; %1' contours /1000000 so they are under the rows
    elseif (max(max(fullsurf(:,:)))-min(min(fullsurf(:,:))))<200
        levels=10000+min(min(fullsurf(:,:)))/1000000:2/1000000:10000+max(max(fullsurf(:,:)))/1000000;
    else
        levels=10000+min(min(fullsurf(:,:)))/1000000:5/1000000:10000+max(max(fullsurf(:,:)))/1000000;
    end
    contour3(surface.xq,surface.yq,10000+fullsurf/1000000,levels,'-k')
    alpha(0.4)
    hold on
    s=surf(surface.xq, surface.yq, foo+10000, colors);
    shading interp
    view(2)
    xlim([min(minx)-50 max(maxx)+50])
    ylim([min(miny)-50 max(maxy)+50])
    s.AlphaData = (foo<=-0.001)|(foo>=0.001);
    s.FaceAlpha = 0.35;
    s.EdgeAlpha = 0.5;
    daspect([1 1 1])
    set(gca, 'units', 'normalized'); %Just making sure it's normalized
    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
    %[Left Bottom Right Top] spacing
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
    set(gca, 'Position', NewPos);
    if const.writefiles == 1
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' cutandfill.png'])
    end

    total_grad.cut=foom;
    total_grad.fill=foop;
    formatSpec = '%.1f';
    tcut = num2str(total_grad.cut,formatSpec);
    tfill = num2str(total_grad.fill,formatSpec);
    tarea = num2str(grd_area,formatSpec);

    import mlreportgen.ppt.*
    t=string(datetime("today"));
    pbt = table(string(tcut),string(tfill),string(tarea),'VariableNames',["Cut (CY)","Fill (CY)","Disturbed Area (Acres)"]);
    ppt = Presentation([const.outpath '/figures/'  char(const.customer) '_' char(const.project)...
        '_GRADING_ESTIMATES_' char(t) '.pptx'],'PPT/grading_fig.potx');
    open(ppt);
    PMN=Picture([const.outpath '/figures/' char(const.sect_name) ' cutandfill.png']);
    replace(ppt,'PlanarMain',PMN)
    %grading brackets
    formatSpec = ' PROJECT IS GRADED TO MAINTAIN REVEAL RANGE OF %.1f FEET TO %.1f FEET ';
    gbstrh=sprintf(formatSpec,const.min_wp,const.max_wp);
    gbstr=Paragraph(gbstrh);
    gbstr.Style ={HAlign('center')};
    gbstr.Font="Hack";
    gbstr.FontSize="8";
    replace(ppt,"GradeBracket",gbstr)
    %project name
    formatSpec = 'PROJECT\n%s';
    pj=upper(string(const.sect_name));
    projnameh=sprintf(formatSpec, pj);
    projname = Paragraph(projnameh);
    projname.Style = {HAlign('center')};
    projname.Font="Hack";
    projname.FontSize="12";
    replace(ppt,"ProjectCoSt",projname)
    %customer
    formatSpec = '%s';
    cn=upper(string(const.customer));
    custnameh=sprintf(formatSpec, cn);
    custname = Paragraph(custnameh);
    custname.Style = {HAlign('center')};
    custname.Font="Hack";
    custname.FontSize="12";
    replace(ppt,"Customer",custname)
    %customer
    formatSpec = '%s\n%s';
    city=upper(string(const.city));
    co_st=upper(string(const.county));
    citycoh=sprintf(formatSpec, city, co_st);
    cityco = Paragraph(citycoh);
    cityco.Style = {HAlign('center')};
    cityco.Font="Hack";
    cityco.FontSize="12";
    replace(ppt,"Coordinates",cityco)


    %revision table
    revtable=Table({'Rev','Date','By','Checked','Approved','Description';...
        const.rev_num, t, "OR", "JR", "OR","AutoGrade";...
        " ", " ", " ", " "," "," "});
    revtable.Font="Hack";
    revtable.FontSize="6";
    revtable.entry(1,1).Style={ColWidth('.35in')};
    revtable.entry(1,2).Style={ColWidth('.7in')};
    revtable.entry(1,3).Style={ColWidth('.30in')};
    revtable.entry(1,4).Style={ColWidth('.5in')};
    revtable.entry(1,5).Style={ColWidth('.55in')};
    revtable.entry(1,6).Style={ColWidth('1.02in')};
    revtable.row(1).Style = [revtable.row(1).Style {RowHeight("0.06in")}];
    revtable.row(2).Height="0.06in";
    revtable.row(3).Height="0.06in";
    replace(ppt,'RevTable',revtable)

    %PileTable
    PileTable=Table(pbt);
    PileTable.Font="Hack";
    PileTable.FontSize="6";
    PileTable.entry(1,1).Style={ColWidth('0.8in')};
    PileTable.entry(1,2).Style={ColWidth('0.8in')};
    PileTable.entry(1,2).Style={ColWidth('0.8in')};
    PileTable.entry(1,2).Style={ColWidth('0.8in')};
    PileTable.row(1).Style = [PileTable.row(1).Style {RowHeight("0.06in")}];
    for j=2:size(pbt,1)
        PileTable.row(j).Height="0.06in";
    end
    replace(ppt,'PileTable',PileTable)
    close(ppt);
    pause(2);

    gname=[const.outpath '/figures/'  char(const.customer) '_' char(const.project) '_GRADING_ESTIMATES_' char(t) '.pptx'];
    delete([const.outpath '/figures/' char(const.sect_name) ' cutandfill.png']);
    if ispc
        pptview(gname,'converttopdf')
        gnamepdf=strrep(gname,'pptx','pdf');
    else
        gnamepdf=gname;
    end
end
end