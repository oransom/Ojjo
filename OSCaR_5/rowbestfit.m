function [kp] = rowbestfit(kp, kpv, F_og, const, surface)
% clc
% clf
% close all

neighbor_distance=const.neighbor_distance;
rpf=const.row_pfit;
npf=const.neigh_pfit;
neighs=const.num_neigh;

if neighs>=1
    neighbor_distance=neighs*neighbor_distance;
else
    fprintf('INFORMATIONAL POPUP: You have chosen no nearest neighbor solution \n')
end

for i=1:length(kpv)
    maxwidth=max(length(kp.(kpv{i}).yi:-const.bin_y:kp.(kpv{i}).yf)); %find max width of the binned y
    kh.(kpv{i}).y_inter=NaN(size(kp.(kpv{i}).yi,1),maxwidth); %create a nan matrix to absorb the whole size
    kh.(kpv{i}).x_inter=NaN(size(kp.(kpv{i}).yi,1),maxwidth); %create a nan matrix to absorb the whole size
    kh.(kpv{i}).z_inter=NaN(size(kp.(kpv{i}).yi,1),maxwidth); %create a nan matrix to absorb the whole size
    for j=1:size(kp.(kpv{i}).yi,1)
        y_int_len=length(kp.(kpv{i}).yi(j):-const.bin_y:kp.(kpv{i}).yf(j)); %find the length of the row
        kh.(kpv{i}).y_inter(j,1:y_int_len)=kp.(kpv{i}).yi(j):-const.bin_y:kp.(kpv{i}).yf(j); %if this errors it might be because y is increasing to the south
        kh.(kpv{i}).x_inter(j,1:y_int_len)=kp.(kpv{i}).xi(j);
    end
    kh.(kpv{i}).z_inter=F_og(kh.(kpv{i}).x_inter, kh.(kpv{i}).y_inter);

    fn=fieldnames(kh.(kpv{i})); %this whole block finds the last vals of each of the row structs (x,y,z) etc
    for fnx=1:numel(fn)
        B=~isnan(kh.(kpv{i}).(fn{fnx}));
        Indices=arrayfun(@(x) find(B(x,:), 1, 'last'), 1:size(kh.(kpv{i}).(fn{fnx}),1));
        lastval.(fn{fnx}) = arrayfun(@(x,y) kh.(kpv{i}).(fn{fnx})(y,x), Indices, 1:size(kh.(kpv{i}).(fn{fnx}), 1));
    end
    
    %polyfit of ns row surface
    for j=1:size(kh.(kpv{i}).y_inter,1) %this is finding the poly=1 bestfit of each of the rows, with elevations interpreted to the bin_y level
        kh.(kpv{i}).c(j,:)=polyfit(kh.(kpv{i}).y_inter(j,1:Indices(j)),kh.(kpv{i}).z_inter(j,1:Indices(j)),1);
        kh.(kpv{i}).f(j,:)=polyval(kh.(kpv{i}).c(j,:),kh.(kpv{i}).y_inter(j,:));
    end
    
    
    kh.(kpv{i}).xnorth=kh.(kpv{i}).x_inter(:,1); %northern value of z
    kh.(kpv{i}).xsouth=arrayfun(@(x,y) kh.(kpv{i}).x_inter(y,x), Indices, 1:size(kh.(kpv{i}).x_inter, 1)); %southern value of z
    kh.(kpv{i}).xsouth=kh.(kpv{i}).xsouth';
    kh.(kpv{i}).ns_comb_x=[kh.(kpv{i}).xnorth; kh.(kpv{i}).xsouth]; %combined n/s for use in fitter
    
    kh.(kpv{i}).ynorth=kh.(kpv{i}).y_inter(:,1); %northern value of z
    kh.(kpv{i}).ysouth=arrayfun(@(x,y) kh.(kpv{i}).y_inter(y,x), Indices, 1:size(kh.(kpv{i}).y_inter, 1)); %southern value of z
    kh.(kpv{i}).ysouth=kh.(kpv{i}).ysouth';
    kh.(kpv{i}).ns_comb_y=[kh.(kpv{i}).ynorth; kh.(kpv{i}).ysouth]; %combined n/s for use in fitter
    
    kh.(kpv{i}).znorth=kh.(kpv{i}).f(:,1); %northern value of z
    kh.(kpv{i}).zsouth=arrayfun(@(x,y) kh.(kpv{i}).f(y,x), Indices, 1:size(kh.(kpv{i}).f, 1)); %southern value of z
    kh.(kpv{i}).zsouth=kh.(kpv{i}).zsouth';
    kh.(kpv{i}).ns_comb_z=[kh.(kpv{i}).znorth; kh.(kpv{i}).zsouth]; %combined n/s for use in fitter

    kh.(kpv{i}).xval=[kh.(kpv{i}).x_inter(:,1); lastval.x_inter']; %this two lines are concat first & last row values
    kh.(kpv{i}).yval=[kh.(kpv{i}).y_inter(:,1); lastval.y_inter'];
    kh.(kpv{i}).xyn=[ones(size(kh.(kpv{i}).x_inter,1),1);2*ones(size(kh.(kpv{i}).x_inter,1),1)]; %assigning 1's for np and 2's for sp
    
    rsx=[kh.(kpv{i}).xval, kh.(kpv{i}).yval];
    rsy=[kh.(kpv{i}).xval, kh.(kpv{i}).yval];
    [Idx,D]=rangesearch(rsx,rsy,neighbor_distance); %search for all row neighbors within neighbor_distance
    maxLengthCell=max(cellfun('size',Idx,2));  %finding the longest vector in the cell array
    for m=1:length(Idx)
        for n=cellfun('size',Idx(m),2)+1:maxLengthCell
            Idx{m}(n)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
        end
    end
    neighbors=cell2mat(Idx); %turning IDX into numerical values
    neighbors(neighbors==0)=NaN; %turning 0s into NaN's so we can avoid working with them

    [m,n]=size(neighbors);
    kh.(kpv{i}).rnbrw=zeros(m,1);
    kh.(kpv{i}).rnbrc=1:1:m;
    kh.(kpv{i}).rnbrc=kh.(kpv{i}).rnbrc';
    kh.(kpv{i}).rnbre=zeros(m,1);
    if neighs==0
        for k=1:m
            for l=2:n
                if ~isnan(neighbors(k,l))
                    if kh.(kpv{i}).xval(neighbors(k,l))<kh.(kpv{i}).xval(neighbors(k,1))...
                            && abs(kh.(kpv{i}).yval(neighbors(k,l))-kh.(kpv{i}).yval(neighbors(k,1)))<3 %finding neigbors to west
                        kh.(kpv{i}).rnbrw(k,1)=neighbors(k,l);
                    elseif kh.(kpv{i}).xval(neighbors(k,l))>kh.(kpv{i}).xval(neighbors(k,1))...
                            && abs(kh.(kpv{i}).yval(neighbors(k,l))-kh.(kpv{i}).yval(neighbors(k,1)))<3 % finding neighbors to east
                        kh.(kpv{i}).rnbre(k,1)=neighbors(k,l);
                    end
                end
            end
        end
    else
        for k=1:m
            for l=2:n
                if ~isnan(neighbors(k,l))
                    if kh.(kpv{i}).xval(neighbors(k,l))<kh.(kpv{i}).xval(neighbors(k,1))...
                            && abs(kh.(kpv{i}).yval(neighbors(k,l))-kh.(kpv{i}).yval(neighbors(k,1)))<3 %finding neigbors to west
                        kh.(kpv{i}).rnbrw(k,l-1)=neighbors(k,l);
                    elseif kh.(kpv{i}).xval(neighbors(k,l))>kh.(kpv{i}).xval(neighbors(k,1))...
                            && abs(kh.(kpv{i}).yval(neighbors(k,l))-kh.(kpv{i}).yval(neighbors(k,1)))<3 % finding neighbors to east
                        kh.(kpv{i}).rnbre(k,l-1)=neighbors(k,l);
                    end
                end
            end
        end
    end
    kh.(kpv{i}).rnbrw(kh.(kpv{i}).rnbrw==0)=NaN;
    kh.(kpv{i}).rnbre(kh.(kpv{i}).rnbre==0)=NaN;

    kh.(kpv{i}).zfin=kh.(kpv{i}).ns_comb_z; 

%     for j=1:size(kh.(kpv{i}).rnbrw,1) %this loop fits all np and sp to each other based on the variable row_fit
%         if ~isnan(kh.(kpv{i}).rnbrw(j)) && ~isnan(kh.(kpv{i}).rnbre(j)) %if there are neighbors to the east and west
%             if kh.(kpv{i}).xyn(kh.(kpv{i}).rnbrw(j)) ==  kh.(kpv{i}).xyn(kh.(kpv{i}).rnbre(j)) %if all north or all south
%                 xf=[kh.(kpv{i}).ns_comb_x(kh.(kpv{i}).rnbrw(j)),kh.(kpv{i}).ns_comb_x(kh.(kpv{i}).rnbre(j))]; %combining into two rows for use in polyfit
%                 zf=[kh.(kpv{i}).ns_comb_z(kh.(kpv{i}).rnbrw(j)),kh.(kpv{i}).ns_comb_z(kh.(kpv{i}).rnbre(j))]; %combining into two rows for use in polyfit
%                 ewpf=polyfit(xf,zf,1); %linear polyfit
%                 zfintemp=polyval(ewpf,kh.(kpv{i}).ns_comb_x(kh.(kpv{i}).rnbrc(j)));
%                 kh.(kpv{i}).zfin(j)=kh.(kpv{i}).ns_comb_z(j)-(row_fit*(kh.(kpv{i}).ns_comb_z(j)-zfintemp)); %check this line: is it plus or minus? 
%             end
%         end
%     end
if neighs>1
    for j=1:size(kh.(kpv{i}).rnbrw,1)
        if neighs>1
            allneighs=[kh.(kpv{i}).rnbrc,~isnan(kh.(kpv{i}).rnbrw).*(kh.(kpv{i}).rnbrw),~isnan(kh.(kpv{i}).rnbre).*(kh.(kpv{i}).rnbre)];
        end
    end

    allneighs=allneighs';
    allneighs(isnan(allneighs))=0; %this block gets rid of extra nans
        T=allneighs~=0;
        n=sum(T,1);
        m=max(n);
        allneighs2=nan(m,size(allneighs,2));
        for k=1:size(allneighs,2)
            allneighs2(1:n(k),k) = allneighs(T(:,k),k);
        end
    allneighs=allneighs2';
    B=~isnan(allneighs);
    Indices=arrayfun(@(x) find(B(x,:), 1, 'last'), 1:size(allneighs,1));
    %lastval.(fn{fnx}) = arrayfun(@(x,y) plf.(fn{fnx})(x,y), Indices, 1:size(plf.(fn{fnx}), 2));
    
    for j=1:size(allneighs,1)
        zpf_neigh=polyfit(kh.(kpv{i}).ns_comb_x(allneighs(j,1:Indices(j))),kh.(kpv{i}).ns_comb_z(allneighs(j,1:Indices(j))),npf); %interpolating along the x/z line of unique y values above
        zpf2_neigh=polyval(zpf_neigh,kh.(kpv{i}).ns_comb_x(allneighs(j,1:Indices(j)))); %creating northern z values along the x northern line of each unique y value 
        for n=1:size(zpf2_neigh)
            zplotn1(j)=zpf2_neigh(1);
        end
    end
end
%the point of this section is to try and establish mean z values for n/s of
%every row
    mask1=find(kh.(kpv{i}).xyn==1);
    mask2=find(kh.(kpv{i}).xyn==2);
    ny=unique(kh.(kpv{i}).ns_comb_y(mask1));
    sy=unique(kh.(kpv{i}).ns_comb_y(mask2));
    zplot3=kh.(kpv{i}).ns_comb_z;
    zplot4=kh.(kpv{i}).ns_comb_z;

    for j=1:(size(ny,1))
        zfoo=find(kh.(kpv{i}).ns_comb_y==ny(j)); %finding unique y values so we can interpolate along that line
        zpf1=polyfit(kh.(kpv{i}).ns_comb_x(zfoo),kh.(kpv{i}).ns_comb_z(zfoo),1); %interpolating along the x/z line of unique y values above
        zpf2=polyval(zpf1,kh.(kpv{i}).ns_comb_x(zfoo)); %creating northern z values along the x northern line of each unique y value 
        zpf1_2=polyfit(kh.(kpv{i}).ns_comb_x(zfoo),kh.(kpv{i}).ns_comb_z(zfoo),rpf); %polyfitting to chosen factor
        zpf2_2=polyval(zpf1_2,kh.(kpv{i}).ns_comb_x(zfoo)); %polyfitting to chosen factor
        for n=1:size(zpf2,1)
            zplot3(zfoo(n))=zpf2(n);
            zplot4(zfoo(n))=zpf2_2(n);
        end
    end
    for j=1:(size(sy,1))
        zfoo=find(kh.(kpv{i}).ns_comb_y==sy(j));
        zpf1=polyfit(kh.(kpv{i}).ns_comb_x(zfoo),kh.(kpv{i}).ns_comb_z(zfoo),1);
        zpf2=polyval(zpf1,kh.(kpv{i}).ns_comb_x(zfoo));
        zpf1_2=polyfit(kh.(kpv{i}).ns_comb_x(zfoo),kh.(kpv{i}).ns_comb_z(zfoo),rpf);
        zpf2_2=polyval(zpf1_2,kh.(kpv{i}).ns_comb_x(zfoo));
        for n=1:size(zpf2,1)
            zplot3(zfoo(n))=zpf2(n);
            zplot4(zfoo(n))=zpf2_2(n);
        end
    end

    for j=1:size(mask1,1)
        kh.(kpv{i}).xplot(j,:)=[kh.(kpv{i}).ns_comb_x(mask1(j)),kh.(kpv{i}).ns_comb_x(mask2(j))];
        kh.(kpv{i}).yplot(j,:)=[kh.(kpv{i}).ns_comb_y(mask1(j)),kh.(kpv{i}).ns_comb_y(mask2(j))];
        kh.(kpv{i}).zplot(j,:)=[kh.(kpv{i}).zfin(mask1(j)),kh.(kpv{i}).zfin(mask2(j))];
        kh.(kpv{i}).zplot3(j,:)=[zplot3(mask1(j)),zplot3(mask2(j))];
        kh.(kpv{i}).zplot4(j,:)=[zplot4(mask1(j)),zplot4(mask2(j))];
        if neighs>1
            kh.(kpv{i}).zplotn(j,:)=[zplotn1(mask1(j)),zplotn1(mask2(j))];
        end
    end

if const.sol_type==0
    levels=min(min(surface.zog(:,:))):1:max(max(surface.zog(:,:))); %1' contours
    figure
    set(gcf, 'Position',  [500, 500, 2400, 400])
    subplot(1,5,1)
    contour3(surface.xq,surface.yq,surface.zog,levels)
    hold on
    for j=1:size(kh.(kpv{i}).y_inter,1)
        plot3(kh.(kpv{i}).x_inter(j,:),...
            kh.(kpv{i}).y_inter(j,:),...
            kh.(kpv{i}).f(j,:),'k-', 'LineWidth',3)
    end
    view(2)
    title('terrain')

    subplot(1,5,2)
    hold on
    for j=1:size(kh.(kpv{i}).y_inter,1)
        plot3(kh.(kpv{i}).x_inter(j,:),...
            kh.(kpv{i}).y_inter(j,:),...
            kh.(kpv{i}).f(j,:),'b-', 'LineWidth',2)
    end
    %zlim([3200 3240])
    view(-51,7)
    title('terrain best fit soluton: Sol 1')

    subplot(1,5,3)
    hold on
    for j=1:size(kh.(kpv{i}).xplot,1)
        plot3(kh.(kpv{i}).xplot(j,:),...
            kh.(kpv{i}).yplot(j,:),...
            kh.(kpv{i}).zplot3(j,:),'r-', 'LineWidth',2)
    end
    %zlim([3200 3240])
    view(-51,7)
    title('linear solution: Sol 2')

    subplot(1,5,4)
    hold on
    for j=1:size(kh.(kpv{i}).xplot,1)
        plot3(kh.(kpv{i}).xplot(j,:),...
            kh.(kpv{i}).yplot(j,:),...
            kh.(kpv{i}).zplot4(j,:),'r-', 'LineWidth',2)
    end
    %zlim([3200 3240])
    view(-51,7)
    title('polynomal whole row solution: Sol 3')
    if neighs>1
        subplot(1,5,5)
        %contour3(surface.xq,surface.yq,surface.zog,levels)
        hold on
        for j=1:size(kh.(kpv{i}).xplot,1)
            plot3(kh.(kpv{i}).xplot(j,:),...
                kh.(kpv{i}).yplot(j,:),...
                kh.(kpv{i}).zplotn(j,:),'r-', 'LineWidth',2)
        end
        %zlim([3200 3240])
        view(-51,7)
        title('neighbor solution: Sol 4')
    end
end
if const.sol_type==0
    prompt = "Choose Solution Type:  ";
    const.sol_type = input(prompt);
end
end
if const.sol_type==1
for j=1:length(kpv)
    kp.(kpv{j}).zi=kh.(kpv{j}).znorth;
    kp.(kpv{j}).zf=kh.(kpv{j}).zsouth;
    kp.(kpv{j}).slp=atand((kp.(kpv{j}).zf-kp.(kpv{j}).zi)./abs(ysouth-ynorth));
end
elseif const.sol_type==2
for j=1:length(kpv)
    kp.(kpv{j}).zi=kh.(kpv{j}).zplot3(:,1);
    kp.(kpv{j}).zf=kh.(kpv{j}).zplot3(:,2);
    kp.(kpv{j}).slp=atand((kp.(kpv{j}).zf-kp.(kpv{j}).zi)./abs(kh.(kpv{j}).ynorth-kh.(kpv{j}).ysouth));
end
elseif const.sol_type==3
for j=1:length(kpv)
    kp.(kpv{j}).zi=kh.(kpv{j}).zplot4(:,1);
    kp.(kpv{j}).zf=kh.(kpv{j}).zplot4(:,2);
    kp.(kpv{j}).slp=atand((kp.(kpv{j}).zf-kp.(kpv{j}).zi)./abs(kh.(kpv{j}).ynorth-kh.(kpv{j}).ysouth));
end
elseif const.sol_type==4
for j=1:length(kpv)
    kp.(kpv{j}).zi=kh.(kpv{j}).zplotn(:,1);
    kp.(kpv{j}).zf=kh.(kpv{j}).zplotn(:,2);
    kp.(kpv{j}).slp=atand((kp.(kpv{j}).zf-kp.(kpv{j}).zi)./abs(kh.(kpv{j}).ynorth-kh.(kpv{j}).ysouth));
end
else
    fprintf('Not a valid choice.   \n')
end

foo=1;
end