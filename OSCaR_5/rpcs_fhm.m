function rpcs_fhm(kpv,pd_out,F_og,surface,ppp_out,const,grd_ext)
%change some variable names for compactness - not best for memory, but
%whatever
p=pd_out; %tpx=zeros(200,200);
fhmpc=readtable("../config/fhm_plot_config.xlsx");
fhmpc.color=([fhmpc.c1/255, fhmpc.c2/255, fhmpc.c3/255]);
dat_bins=const.cott_bins{1};
%create RGB's for dat_bins
for j=1:length(dat_bins.bin_color)
    if matches(string(dat_bins.bin_color(j)),"b")
        dat_bins.rgb{j}=[0 0 1];
    elseif matches(string(dat_bins.bin_color(j)),"g")
        dat_bins.rgb{j}=[0 1 0];
    elseif matches(string(dat_bins.bin_color(j)),"y")
        dat_bins.rgb{j}=[1 1 0];
    elseif matches(string(dat_bins.bin_color(j)),"o")
        dat_bins.rgb{j}=[1 .4 0];
    else
        dat_bins.rgb{j}=[1 0 0];
    end
end

yi=200;
xi=21; %should really be an odd number
xpdim=const.xpdim; %panel height in mm
xm(:)=(-xpdim/25.4/2:xpdim/25.4/(xi-1):xpdim/25.4/2)/12;
offset=0;%const.offset; %.33 for ati, 2/12 for nevados 1.5/12 NXT
max_tilt=const.maxtilt; %degrees 52 ati 60 nxt
leadedge=const.leading_edge; %leading edge clearance

%rollup some data
for i=1:length(kpv)
    if i==1
        cott_ft_ru=pd_out.(kpv{i}).cott_ft;
        pile_reveal_ru=pd_out.(kpv{i}).pile_reveal;
        slp_deg_ru=pd_out.(kpv{i}).slp_deg;
    else
        cott_ft_ru=[cott_ft_ru;pd_out.(kpv{i}).cott_ft];
        pile_reveal_ru=[pile_reveal_ru;pd_out.(kpv{i}).pile_reveal];
        slp_deg_ru=[slp_deg_ru;pd_out.(kpv{i}).slp_deg];
    end
end

figure ('Visible', 'off')
set(gcf, 'Position',  [100, 100, 2400, 1800])
levels=min(min(surface.zog(:,:)))/1000000:1/1000000:max(max(surface.zog(:,:)))/1000000; %1' contours /1000000 so they are under the rows
contour3(surface.xq,surface.yq,surface.zog/1000000,levels,':k');
alpha(0.3)
hold on
view(2)
%% ANGLE == FLAT
for i=1:length(kpv)
    numsect=1:max(p.(kpv{i}).section);
    for j=1:max(numsect)
        numrow(j)=max(p.(kpv{i}).row_number(p.(kpv{i}).section==j));
    end
    l=0;
    for j=1:max(numsect)
        for k=1:max(numrow(j))
            l=l+1;
            tpx{:,l}=p.(kpv{i}).tpx(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
            tpy{:,l}=p.(kpv{i}).tpy(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
            top{:,l}=p.(kpv{i}).cott_ft(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
            elev{:,l}=p.(kpv{i}).tpz(p.(kpv{i}).section==j & p.(kpv{i}).row_number==k);
        end
    end
    minx(i)=min(p.(kpv{i}).tpx(:));
    maxx(i)=max(p.(kpv{i}).tpx(:));
    miny(i)=min(p.(kpv{i}).tpy(:));
    maxy(i)=max(p.(kpv{i}).tpy(:));

    %if a error happens somewhere around here check that each section is restarting at 1 and not continuing
    for j=1:length(tpx)
        if ~isempty(tpx{j})
            x=cell2mat(tpx(j));
            y=cell2mat(tpy(j));
            h=cell2mat(top(j));
            e=cell2mat(elev(j));
            y_int=repmat(y(1):(y(end)-y(1))/yi:y(end),xi,1);
            yh(:,:,j)=y_int;
            x_int=ones(size(y_int)).*(xm+mean(x))';
            xh(:,:,j)=x_int;
            e_int=interp1(y,e,y_int);
            eh(:,:,j)=e_int;
        end
    end
    for j=1:length(tpx)
        if ~isempty(tpx{j})
            hag=eh(:,:,j)-F_og(xh(:,:,j),yh(:,:,j)); %for future FYI we can use histcounts here with edges to bin distance above ground
            panels.(kpv{i}).hagh(:,:,j)=hag+offset;
        end
    end

    for j=1:length(tpx)
        if ~isempty(tpx{j})
            Zhold=reshape(panels.(kpv{i}).hagh(:,:,j),[],1);
            X=reshape(xh(:,:,j),[],1);
            Y=reshape(yh(:,:,j),[],1);
            for k=1:length(dat_bins.bin_start)
                Z=zeros(length(Zhold),1);
                if k==1
                    Z_ind=panels.(kpv{i}).hagh(:,:,j)<=dat_bins.bin_finish(k);
                elseif k==length(dat_bins.bin_start)
                    Z_ind=panels.(kpv{i}).hagh(:,:,j)>dat_bins.bin_start(k);
                else
                    Z_ind=panels.(kpv{i}).hagh(:,:,j)>dat_bins.bin_start(k) & panels.(kpv{i}).hagh(:,:,j)<=dat_bins.bin_finish(k);
                end
                Z(Z_ind)=Zhold(Z_ind);
                Z(Z==0)=NaN;
                s1=scatter3(X,Y,Z,5,dat_bins.rgb{k});
                alpha(s1,.25)
                clear Z
            end
        end
    end
    for j=1:length(fhmpc.c1)
        fp{j}=find(contains(p.(kpv{i}).pt_desc,string(fhmpc.type(j))) & contains(p.(kpv{i}).xsect,string(fhmpc.xsect(j))));
        if ~isempty(fp{j})
            scatter3(p.(kpv{i}).tpx(fp{j}),p.(kpv{i}).tpy(fp{j}),p.(kpv{i}).tpz(fp{j}),700,...
                [fhmpc.c1(j)/255,fhmpc.c2(j)/255,fhmpc.c3(j)/255],"filled",string(fhmpc.symbol(j)))
            text(p.(kpv{i}).tpx(fp{j})-2.75,p.(kpv{i}).tpy(fp{j}),p.(kpv{i}).tpz(fp{j})+1,string(p.(kpv{i}).total_steel(fp{j})))
        end
    end


    %GRADING EXTENT PLOTTER
    if ~isempty(grd_ext)
    for j=1:length(grd_ext.gbx)
        plot(grd_ext.xgh{j}(grd_ext.gbx{j}),grd_ext.ygh{j}(grd_ext.gby{j}),'k--','LineWidth',2)
        hold on
    end
    end

    axis off
    hold on
    xlim([min(minx)-50 max(maxx)+50])
    ylim([min(miny)-50 max(maxy)+50])

    set(gca, 'units', 'normalized'); %Just making sure it's normalized
    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
    %[Left Bottom Right Top] spacing
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
    set(gca, 'Position', NewPos);
    if const.min_wp ~= const.max_wp
        clim([const.min_wp const.max_wp])
    else
        clim([const.min_wp+offset-.5 const.max_wp+offset+.5])
    end
end
if const.writefiles == 1
    if const.simple_out==0
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' fhm_plot.png'])
    else
        mkdir([const.outpath '/figures' ]);
        saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' fhm_plot.png'])
    end
    foo=1;
end


