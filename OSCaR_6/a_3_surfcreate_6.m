function [surface] = a_3_surfcreate_6(const)

gtopo=readtable(string(const.topo)); %read in topo
gtopo=rmmissing(gtopo); %remove any missing or bad topo values
surface.F_og=scatteredInterpolant(gtopo.Var1,gtopo.Var2,gtopo.Var3,'natural','none'); %queryable surface, natural internal interp, no extrapolation
[surface.xq,surface.yq]=meshgrid(min(gtopo.Var1):const.bin_x:max(gtopo.Var1), min(gtopo.Var2):const.bin_y:max(gtopo.Var2));
surface.zog=surface.F_og(surface.xq,surface.yq);

if ~strcmpi(const.flood,'na') & ~isnumeric(const.flood) %if there is a flood layer and the flood layer is not a flood depth
    flood=readtable(string(const.flood));  %read in the flood layer
    surface.Flood_s=scatteredInterpolant(flood.Var1,flood.Var2,flood.Var3,'natural','none'); %create the flood surface
    surface.flood=surface.Flood_s(surface.xq,surface.yq); %sample the flood surface on the same grid as the topo
    surface.flood(surface.flood<surface.zog+0.01)=NaN;%surface.zog(surface.flood<surface.zog+0.01); %any flood element that is below surface element becomes surface element
    if isequal(surface.flood,surface.zog) %if entirety of flood surface is below surface, we assume it's a flood depth layer
        surface.flood(surface.flood<0.01)=0; %close to zero cells become zero
        surface.flood=surface.zog+surface.flood; %flood surface is added to surface
    end
    surface.floodmask=isnan(surface.flood); %create floodmask
    surface.Flood_s=scatteredInterpolant(reshape(surface.xq,[],1),...
        reshape(surface.yq,[],1),reshape(surface.flood,[],1),'natural','none'); %recreate new queryable Flood surface
elseif ~strcmpi(const.flood,'na') & isnumeric(const.flood) %if a numeric value is given set everything below that value to flood depth given
    flooded=surface.zog<const.flood; %find surface below given flood depth
    dry=surface.zog>const.flood; %find surface above flood depth
    surface.flood=zeros(size(surface.zog));
    surface.flood(flooded)=const.flood;
    surface.flood(dry)=surface.zog(dry);
    surface.Flood_s=scatteredInterpolant(reshape(surface.xq,[],1),...
        reshape(surface.yq,[],1),reshape(surface.flood,[],1),'natural','none'); %recreate new queryable Flood surface
    surface.floodmask=dry; %create floodmask
else
    surface.Flood_s=[];
end

end

