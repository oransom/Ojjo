function [cpl, uc, trs, pmd, oos_data, ns_ugly] = pilelocationcalc(F_og,Flood_s,kpv,dat,span,const,surface)
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

if const.slope_or_bf==1%next five lines were just added
    [kp] = slopefit(kp, kpv, F_og, const); %I need to figure out how to move this out of here
elseif const.slope_or_bf==2
    [kp] = rowbestfit(kp, kpv, F_og, const, surface);
end

for i=1:length(kpv)
    if const.slope_or_bf==2
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
    pmd.(kpv{i}).ismotor=contains(cellstr(pmd.(kpv{i}).tt),'mtr');

    [~,n]=size(kp.(kpv{i}).span);
    pint=[kp.(kpv{i}).xi,kp.(kpv{i}).yi,cpl.(kpv{i}).zit]; %pile initial top
    if const.slope_or_bf==1
        pf=[kp.(kpv{i}).xf,kp.(kpv{i}).yf,cpl.(kpv{i}).zft];  %pile final top
    else
        pf=[kp.(kpv{i}).xf,kp.(kpv{i}).yfa,cpl.(kpv{i}).zft];  %pile final top
    end
    uv=(pf-pint)./vecnorm(pf-pint,2,2); %find unit vector of all rows along each row of pf-pint (2 norm of dimension 2 [rows])
    cpl.(kpv{i}).tpx=pint(:,1)-kp.(kpv{i}).span.*uv(:,1);
    cpl.(kpv{i}).tpy=pint(:,2)-kp.(kpv{i}).span.*uv(:,2);
    cpl.(kpv{i}).tpz=pint(:,3)-kp.(kpv{i}).span.*uv(:,3);
    cpl.(kpv{i}).slp=repmat(kp.(kpv{i}).slp,1,n);
    pmd.(kpv{i}).span_type=repmat(dat.(kpv{i}).span,1,n);
    cpl.(kpv{i}).bpz=F_og(cpl.(kpv{i}).tpx,cpl.(kpv{i}).tpy); %find z plumb to earth from top of cpl
    cpl.(kpv{i}).minlenf=repmat(kp.(kpv{i}).minlen',1,n);
    cpl.(kpv{i}).mh=cpl.(kpv{i}).tpz-cpl.(kpv{i}).bpz; %find length of all trusses

    if const.slope_or_bf==1
        cpl.(kpv{i}).paddtl=(min(cpl.(kpv{i}).mh,[],2)<kp.(kpv{i}).minlen').*(kp.(kpv{i}).minlen'-min(cpl.(kpv{i}).mh,[],2));
        %pile additional - min in mh row less than minlen for that row then add difference between it and minlegth to tpz to get tpzf
    elseif const.slope_or_bf==2
        if const.reveal_bias==0
            centerline=const.min_wp;
        else
            centerline=const.min_wp+(const.max_wp-const.min_wp)/(1/const.reveal_bias); %i've thought a lot about this and it makes sense
        end
        cpl.(kpv{i}).paddtl=centerline-mean(cpl.(kpv{i}).mh,2,'omitnan'); %brings mh down to centerline
    end

    cpl.(kpv{i}).tpzf=cpl.(kpv{i}).tpz+cpl.(kpv{i}).paddtl; %adding what's needed back to tpz original (TPZF - where the F is for final)
    cpl.(kpv{i}).mhc=cpl.(kpv{i}).tpzf-cpl.(kpv{i}).bpz; %min height check - all cpl heights should be equal/greater than min_wp

    %%THIS WHOLE SECTION IS FLOOD
    %bring up for flood elevation
    if strcmpi(const.flood,'na')==0
        cpl.(kpv{i}).freeboard=cpl.(kpv{i}).mhc-(Flood_s(cpl.(kpv{i}).tpx,cpl.(kpv{i}).tpy)-F_og(cpl.(kpv{i}).tpx,cpl.(kpv{i}).tpy));
        cpl.(kpv{i}).freeboard(isnan(cpl.(kpv{i}).freeboard))=const.freeboard;
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

    %motor alignment
    if const.motorc==1
    mtr_slp_shft=max(abs(cpl.(kpv{i}).bpy.*(pmd.(kpv{i}).ismotor)-kp.(kpv{i}).nsy.*(pmd.(kpv{i}).ismotor)),[],2);
    cpl.(kpv{i}).bpy=cpl.(kpv{i}).bpy-mtr_slp_shft;
    cpl.(kpv{i}).tpy=cpl.(kpv{i}).bpy;
    foo=1;
    end

    %slope limiter
    if const.slopelimit==1
        [ns_lim_out] = nslimiter(cpl.(kpv{i}),const,kp.(kpv{i}),F_og);
        cpl.(kpv{i}).tpx=ns_lim_out.tpx;
        cpl.(kpv{i}).tpy=ns_lim_out.tpy;
        cpl.(kpv{i}).tpzf=ns_lim_out.tpzf;
        cpl.(kpv{i}).slp=ns_lim_out.slp;
        cpl.(kpv{i}).bpz=ns_lim_out.bpz;
        cpl.(kpv{i}).bpy=ns_lim_out.bpy;
        cpl.(kpv{i}).bpx=ns_lim_out.bpx;
        cpl.(kpv{i}).mh=ns_lim_out.mh;
    end

    %plugging this here now to make the function run
    if const.nsfit==1
        [cplout] = nsfitter(cpl.(kpv{i}),const,kp.(kpv{i}),F_og);
        cpl.(kpv{i}).tpx=cplout.tpx;
        cpl.(kpv{i}).tpy=cplout.tpy;
        cpl.(kpv{i}).tpzf=cplout.tpzf;
        cpl.(kpv{i}).slp=cplout.slp;
        cpl.(kpv{i}).bpz=cplout.bpz;
        cpl.(kpv{i}).bpy=cplout.bpy;
        cpl.(kpv{i}).bpx=cplout.bpx;
        cpl.(kpv{i}).mh=cplout.mh;
    end

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

for i=1:length(kpv)
    %cpl.(kpv{i}).tbyd=cpl.(kpv{i}).tpy-cpl.(kpv{i}).bpy;
    %cpl.(kpv{i}).tbxd=cpl.(kpv{i}).tpx-cpl.(kpv{i}).bpx;
    cpl.(kpv{i}).tbyd=cpl.(kpv{i}).bpy-cpl.(kpv{i}).nsy;
    cpl.(kpv{i}).tbxd=cpl.(kpv{i}).bpx-cpl.(kpv{i}).nsx;
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