function [cpl, uc, trs, pmd, oos_data, ns_ugly] = pilelocationcalc_Nev(F_og,Flood_s,kpv,dat,span,const,surface)
% cpl (calculated pile location) table contains all top of cpl and bottom of cpl information
% trs table includes all truss location information
% pile metadata (tt = truss type, tl = truss letter, sect = section, row = row)

for i=1:length(kpv)
    kp.(kpv{i}).xi=dat.(kpv{i}).tpx;
    kp.(kpv{i}).yi=dat.(kpv{i}).tpy;
    kp.(kpv{i}).xf=dat.(kpv{i}).tpx;
    for j=1:length(dat.(kpv{i}).tpx)
        kp.(kpv{i}).yf(j,1)=dat.(kpv{i}).tpy(j)...
            -span.(string(dat.(kpv{i}).span(j))).kps(length(span.(string(dat.(kpv{i}).span(j))).kps));
        kp.(kpv{i}).len(j,1)=span.(string(dat.(kpv{i}).span(j))).kps(length(span.(string(dat.(kpv{i}).span(j))).kps));
    end
end
for i=1:length(kpv)
    if const.slope_or_bf==1
        [kp] = slopefit(kp, kpv, F_og, const);
    elseif const.slope_or_bf==2
        [kp] = rowbestfit(kp, kpv, F_og, const, surface);
        kp.(kpv{i}).minlen(kp.(kpv{i}).slp>=const.wp_slope_range | kp.(kpv{i}).slp<=-const.wp_slope_range)=const.wp_slope_outrange;
        kp.(kpv{i}).minlen(kp.(kpv{i}).slp<const.wp_slope_range & kp.(kpv{i}).slp>-const.wp_slope_range)=const.wp_slope_inrange;
        kp.(kpv{i}).yfa=kp.(kpv{i}).yf+tand(kp.(kpv{i}).slp).*(kp.(kpv{i}).zi-kp.(kpv{i}).zf);
    end

    cpl.(kpv{i}).zit=kp.(kpv{i}).zi+kp.(kpv{i}).minlen'; %new top of initial laser cpl
    cpl.(kpv{i}).zft=kp.(kpv{i}).zf+kp.(kpv{i}).minlen'; %new top of final laser cpl
    for j=1:length(dat.(kpv{i}).tpx)
        for k=1:size(span.(string(dat.(kpv{i}).span(j))).tt,1) %there may be some weird loop sizing thing here to check, but I think it's solid
            kp.(kpv{i}).span(j,k)=-span.(string(dat.(kpv{i}).span(j))).kps(k);
            pmd.(kpv{i}).tt{j,k}=span.(string(dat.(kpv{i}).span(j))).tt{k}; %pile meta data - truss type
            pmd.(kpv{i}).tl{j,k}=span.(string(dat.(kpv{i}).span(j))).tl{k}; %pile meta data - truss letter
        end
    end
    kp.(kpv{i}).span(kp.(kpv{i}).span==0)=NaN;
    kp.(kpv{i}).span(:,1)=0;
    kp.(kpv{i}).nsy=kp.(kpv{i}).yi+kp.(kpv{i}).span; %pile location, no slope, used to check how much cpl are moving due to code
    kp.(kpv{i}).nsx=kp.(kpv{i}).xi.*(kp.(kpv{i}).nsy>0); %pile location, no slope, used to check how much cpl are moving due to code
    kp.(kpv{i}).nsx(kp.(kpv{i}).nsx==0)=NaN;

    pmd.(kpv{i}).sect=dat.(kpv{i}).sect.*(kp.(kpv{i}).span<=0);
    pmd.(kpv{i}).row=dat.(kpv{i}).row.*(kp.(kpv{i}).span<=0);
    %inserting third leg information here
    %check for motor location
    idx = cellfun(@isempty,pmd.(kpv{i}).tt); % Find the indexes of empty cell
    pmd.(kpv{i}).tt(idx) = {'NA'}; %insert NA's in empty cells
    pmd.(kpv{i}).ismotor=~cellfun('isempty',strfind(cellstr(pmd.(kpv{i}).tt),'Motor'));

    spanfoo=kp.(kpv{i}).span; %copying to spanholder
    spanfoo(:,1)=1; %setting first column to 1, so doesn't trip NaN
    zfl=spanfoo';
    zfl(isnan(zfl))=0; %this block gets rid of extra nans
    T=zfl~=0;
    last_s=sum(T,1); %this seems to be the variable we need?
    fprintf('starting iteration loop');  
    [m,n]=size(kp.(kpv{i}).span);
    for j=1:m
        if mod(j,50)==0
            fprintf('you are on nevrow iter %5.0f of %5.0f \n',j,m);
        end
        for k=1:n
            if k==1
                cpl.(kpv{i}).tpx(j,k)=kp.(kpv{i}).xi(j);
                cpl.(kpv{i}).tpy(j,k)=kp.(kpv{i}).yi(j);
                cpl.(kpv{i}).tpz(j,k)=F_og(cpl.(kpv{i}).tpx(j,k),cpl.(kpv{i}).tpy(j,k))+kp.(kpv{i}).minlen(j);
            else
                delz=1;
                cpl.(kpv{i}).tpx(j,k)=kp.(kpv{i}).nsx(j,k);
                zup=F_og(cpl.(kpv{i}).tpx(j,k-1),cpl.(kpv{i}).tpy(j,k-1));
                zdown=F_og(cpl.(kpv{i}).tpx(j,k),kp.(kpv{i}).nsy(j,k));
                yup=cpl.(kpv{i}).tpy(j,k-1);
                ydown=cpl.(kpv{i}).tpy(j,k-1)+(kp.(kpv{i}).span(j,k)-kp.(kpv{i}).span(j,k-1));
                y1=ydown-yup;

                itrct=0;
                while delz>0.01 %set really small for testing
                    itrct=itrct+1;
                    zdel=zup-zdown;
                    if itrct>50
                        fprintf('iteration going to shit delz = %2.5f, row number %5.0f \n',delz,j);
                        break
                    end
                    zdown_old=zdown;
                    ydel=ydown-yup;
                    newspan=-sqrt(abs(ydel^2-zdel^2)); %negative to keep w/ span convention
                    zdown=F_og(cpl.(kpv{i}).tpx(j,k),cpl.(kpv{i}).tpy(j,k-1)+newspan);
                    delz=abs(zdown_old-zdown);
                end
                cpl.(kpv{i}).tpy(j,k)=cpl.(kpv{i}).tpy(j,k-1)+newspan;
                cpl.(kpv{i}).tpz(j,k)=F_og(cpl.(kpv{i}).tpx(j,k),cpl.(kpv{i}).tpy(j,k))+kp.(kpv{i}).minlen(j);
            end
        end
    end

    %find slopes upstream and downstream - south slope is negative, north
    %slope is positive
    for j=1:m
        for k=1:last_s(j)
            if k==1
                cpl.(kpv{i}).sloped(j,k)=atand((cpl.(kpv{i}).tpz(j,k+1)-cpl.(kpv{i}).tpz(j,k))...
                    /abs(cpl.(kpv{i}).tpy(j,k)-cpl.(kpv{i}).tpy(j,k+1)));
                cpl.(kpv{i}).slopeu(j,k)=0;
                cpl.(kpv{i}).sumslope(j,k)=cpl.(kpv{i}).sloped(j,k)-cpl.(kpv{i}).slopeu(j,k);
                %above line for slope is minus because that makes sense w/
                %passthrough slope
            elseif k==last_s(j)
                cpl.(kpv{i}).slopeu(j,k)=atand((cpl.(kpv{i}).tpz(j,k)-cpl.(kpv{i}).tpz(j,k-1))...
                    /abs(cpl.(kpv{i}).tpy(j,k-1)-cpl.(kpv{i}).tpy(j,k)));
                cpl.(kpv{i}).sloped(j,k)=0;
                cpl.(kpv{i}).sumslope(j,k)=cpl.(kpv{i}).sloped(j,k)-cpl.(kpv{i}).slopeu(j,k);
            else
                cpl.(kpv{i}).sloped(j,k)=atand((cpl.(kpv{i}).tpz(j,k+1)-cpl.(kpv{i}).tpz(j,k))...
                    /abs(cpl.(kpv{i}).tpy(j,k)-cpl.(kpv{i}).tpy(j,k+1)));
                cpl.(kpv{i}).slopeu(j,k)=atand((cpl.(kpv{i}).tpz(j,k)-cpl.(kpv{i}).tpz(j,k-1))...
                    /abs(cpl.(kpv{i}).tpy(j,k-1)-cpl.(kpv{i}).tpy(j,k)));
                cpl.(kpv{i}).sumslope(j,k)=cpl.(kpv{i}).sloped(j,k)-cpl.(kpv{i}).slopeu(j,k);
            end
        end
    end
    cpl.(kpv{i}).NevJoint=ones(size(cpl.(kpv{i}).sumslope));
    cpl.(kpv{i}).NevJoint(cpl.(kpv{i}).sumslope<=2.0)=1;
    cpl.(kpv{i}).NevJoint(cpl.(kpv{i}).sumslope>2.0 & cpl.(kpv{i}).sumslope<=7.41)=2;
    cpl.(kpv{i}).NevJoint(cpl.(kpv{i}).sumslope>7.41 & cpl.(kpv{i}).sumslope<=14.57)=3;
    cpl.(kpv{i}).NevJoint(cpl.(kpv{i}).sumslope>14.57)=999;
    cpl.(kpv{i}).slp=repmat(kp.(kpv{i}).slp,1,n);
    pmd.(kpv{i}).span_type=repmat(dat.(kpv{i}).span,1,n);
    cpl.(kpv{i}).bpz=F_og(cpl.(kpv{i}).tpx,cpl.(kpv{i}).tpy); %find z plumb to earth from top of cpl
    cpl.(kpv{i}).minlenf=repmat(kp.(kpv{i}).minlen',1,n);
    cpl.(kpv{i}).mh=cpl.(kpv{i}).tpz-cpl.(kpv{i}).bpz; %find length of all trusses

    if const.slope_or_bf==1
        cpl.(kpv{i}).paddtl=(min(cpl.(kpv{i}).mh,[],2)<kp.(kpv{i}).minlen').*0; %no extra pile for Nevados
    elseif const.slope_or_bf==2
        if const.reveal_bias==0
            centerline=0; %no movement for nevados
        else
            centerline=0;
        end
        cpl.(kpv{i}).paddtl=centerline-mean(cpl.(kpv{i}).mh,2,'omitnan'); %brings mh down to centerline
    end

    cpl.(kpv{i}).tpzf=cpl.(kpv{i}).tpz;%+cpl.(kpv{i}).paddtl;

    cpl.(kpv{i}).mhc=cpl.(kpv{i}).tpzf-cpl.(kpv{i}).bpz; %min height check - all cpl heights should be equal/greater than min_wp

    %%THIS WHOLE SECTION IS FLOOD
    %bring up for flood elevation
    if strcmpi(const.flood,'na')==0
        cpl.(kpv{i}).freeboard=cpl.(kpv{i}).mhc-Flood_s(cpl.(kpv{i}).tpx,cpl.(kpv{i}).tpy);
        cpl.(kpv{i}).paddtl_flood=(min(cpl.(kpv{i}).freeboard,[],2)<const.freeboard).*(const.freeboard-min(cpl.(kpv{i}).freeboard,[],2));
        cpl.(kpv{i}).tpzf=cpl.(kpv{i}).tpzf+cpl.(kpv{i}).paddtl_flood;
    end

    %% Back to regular code
    cpl.(kpv{i}).span=kp.(kpv{i}).span;
    cpl.(kpv{i}).bpx=cpl.(kpv{i}).tpx;
    %uc
    uc.(kpv{i}).bpx=cpl.(kpv{i}).tpx;
    uc.(kpv{i}).bpy=cpl.(kpv{i}).tpy;
    uc.(kpv{i}).tpz=cpl.(kpv{i}).tpz;
    uc.(kpv{i}).tpzf=cpl.(kpv{i}).tpzf;
    uc.(kpv{i}).mhc=cpl.(kpv{i}).mhc;

    if const.perp==1
        cpl.(kpv{i}).bpy=cpl.(kpv{i}).tpy+sind(cpl.(kpv{i}).slp).*cpl.(kpv{i}).mhc; %sind because mhc is hypot
    else
        cpl.(kpv{i}).bpy=cpl.(kpv{i}).tpy;
    end
    cpl.(kpv{i}).bpzf=F_og(cpl.(kpv{i}).bpx,cpl.(kpv{i}).bpy);

    %reshape to columns
    cpl.(kpv{i}).ismotor=pmd.(kpv{i}).ismotor'; %make a cpl entry that is whether or not something is a motor.
    cpl.(kpv{i}).tpzf=cpl.(kpv{i}).tpzf';
    cpl.(kpv{i}).tpx=cpl.(kpv{i}).tpx';
    cpl.(kpv{i}).tpy=cpl.(kpv{i}).tpy';
    cpl.(kpv{i}).mhc=cpl.(kpv{i}).mhc';
    cpl.(kpv{i}).bpx=cpl.(kpv{i}).tpx;
    cpl.(kpv{i}).bpz=cpl.(kpv{i}).bpz';
    cpl.(kpv{i}).bpy=cpl.(kpv{i}).bpy';
    cpl.(kpv{i}).bpzf=cpl.(kpv{i}).bpzf';
    cpl.(kpv{i}).nsy=kp.(kpv{i}).nsy';
    cpl.(kpv{i}).nsx=kp.(kpv{i}).nsx';

    pmd.(kpv{i}).tt=reshape(pmd.(kpv{i}).tt',[],1);
    pmd.(kpv{i}).tl=reshape(pmd.(kpv{i}).tl',[],1);
    pmd.(kpv{i}).sect=reshape(pmd.(kpv{i}).sect',[],1);
    pmd.(kpv{i}).row=reshape(pmd.(kpv{i}).row',[],1);
    pmd.(kpv{i}).ismotor=reshape(pmd.(kpv{i}).ismotor',[],1);
    pmd.(kpv{i}).span_type=reshape(pmd.(kpv{i}).span_type',[],1);
    cpl.(kpv{i}).minlenf=reshape(cpl.(kpv{i}).minlenf,[],1);

    oos_data.(kpv{i}).row_tpx_avg=mean(cpl.(kpv{i}).tpx,"omitnan")';
    oos_data.(kpv{i}).row_tpy_avg=cpl.(kpv{i}).tpy(1,:)'-max(abs(cpl.(kpv{i}).span),[],2,"omitmissing")/2;
    oos_data.(kpv{i}).row_tpz_avg=mean(cpl.(kpv{i}).tpzf,"omitnan")';
    oos_data.(kpv{i}).preveal_mean=mean(cpl.(kpv{i}).mhc,"omitnan")';
    oos_data.(kpv{i}).preveal_max=max(cpl.(kpv{i}).mhc)';
    oos_data.(kpv{i}).preveal_min=min(cpl.(kpv{i}).mhc)';
    oos_data.(kpv{i}).preveal_remain=(const.max_wp-max(cpl.(kpv{i}).mhc))';
    oos_data.(kpv{i}).row_length=max(abs(cpl.(kpv{i}).span),[],2,"omitmissing");
end

%get nevados table corners
nevtablecorns(cpl, kpv, pmd, const);


for i=1:length(kpv)
    cpl.(kpv{i}).tbyd=cpl.(kpv{i}).tpy-cpl.(kpv{i}).bpy;
    cpl.(kpv{i}).tbxd=cpl.(kpv{i}).tpx-cpl.(kpv{i}).bpx;
end

for i=1:length(kpv)
    trs.(kpv{i}).wy=cpl.(kpv{i}).bpy;
    trs.(kpv{i}).ey=cpl.(kpv{i}).bpy;
    if const.thirdleg==1
        trs.(kpv{i}).mx=zeros(size(cpl.(kpv{i}).bpx));
        trs.(kpv{i}).my=zeros(size(cpl.(kpv{i}).bpx));
        trs.(kpv{i}).mz=zeros(size(cpl.(kpv{i}).bpx));
    end
    [m,n]=size(cpl.(kpv{i}).span);
    for j=1:n
        for k=1:m %This is where you you add an if statement to make up for the wider motor mounts, you would also have to change pile_truss_output.m to recognize the wider mount
            trs.(kpv{i}).wx(j,k)=cpl.(kpv{i}).tpx(j,k)-tand(const.truss_leg_angle)*cpl.(kpv{i}).mhc(j,k);
            trs.(kpv{i}).ex(j,k)=cpl.(kpv{i}).tpx(j,k)+tand(const.truss_leg_angle)*cpl.(kpv{i}).mhc(j,k);
            trs.(kpv{i}).wz(j,k)=F_og(trs.(kpv{i}).wx(j,k),trs.(kpv{i}).wy(j,k));
            trs.(kpv{i}).ez(j,k)=F_og(trs.(kpv{i}).ex(j,k),trs.(kpv{i}).ey(j,k));
            if const.thirdleg==1
                if cpl.(kpv{i}).ismotor(j,k)==1
                    syms ys zs
                    %ms is slope of line for zs=ms(ys)+bs
                    %ms1 is a projected line in space to approximate the
                    %local slope - tand(40) could be replaced with any
                    %reasonable guess for the distance the TL is from the
                    %truss
                    ms1=(cpl.(kpv{i}).bpzf(j,k)-F_og(cpl.(kpv{i}).bpx(j,k),cpl.(kpv{i}).bpy(j,k)+tand(40)*cpl.(kpv{i}).mhc(j,k)))/(tand(40)*cpl.(kpv{i}).mhc(j,k));
                    bs1=cpl.(kpv{i}).bpzf(j,k);
                    eqn1 = zs == ms1*ys+bs1;
                    ms2=-tand(90-40+cpl.(kpv{i}).slp(k,1));
                    bs2=cpl.(kpv{i}).tpzf(j,k);
                    eqn2 = zs == ms2*ys+bs2;
                    mt_slope_sol=solve([eqn1, eqn2], [ys, zs]);
                    trs.(kpv{i}).mx(j,k)=cpl.(kpv{i}).bpx(j,k);
                    trs.(kpv{i}).my(j,k)=cpl.(kpv{i}).bpy(j,k)+double(mt_slope_sol.ys);
                    trs.(kpv{i}).mz(j,k)=double(mt_slope_sol.zs);
                end
            end
        end
    end
end

for i=1:length(kpv)
    cpl.(kpv{i}).piletype=zeros(size(cpl.(kpv{i}).tpx));
    B = ~isnan(cpl.(kpv{i}).tpx);
    Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(cpl.(kpv{i}).tpx, 2));
    cpl.(kpv{i}).piletype(1,:)=1;
    [~,n]=size(cpl.(kpv{i}).tpx);
    for j=1:n
        cpl.(kpv{i}).piletype(Indices(1,j),j)=2;
    end
end
%% search for non-pretty areas at N/S end rows
[ns_ugly] = northsouth(kpv, cpl, const);
end