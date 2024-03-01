function [const_all] = a_2_constants_config_6(config_file)
baseFileName = string(config_file);        % using a fully-qualified filename here would be good practice
baseFileName = strrep(baseFileName,'\','/');
unix_cf_path = strrep(baseFileName,'\','/');
str=unix_cf_path;
name = wildcardPattern("Except","/");
pat = "/" + name + textBoundary;
outpath = extractBefore(str,pat);
sheetNames = sheetnames(baseFileName);
sheetNames(1) = [];
for k = 1 : numel(sheetNames)
    const_sheet{k} = readcell(baseFileName,'Sheet',sheetNames{k}); %read in each sheet individually
end
for k=1:numel(const_sheet)
    if  k<=3
        cs_trans=const_sheet{k}'; %transpose each sheet
        writecell(cs_trans,'foo.xlsx'); %write the transposed single sheet table back to an excel file
        ct=readtable('foo.xlsx'); %read that same transposed table back in
        delete foo.xlsx %delete the holder excel file
        ct(:,1)=[]; %remove original column headers from variable names
    end
    const_shts.(sheetNames(k))=ct; %create a structure of sheets with everything we need for the job
end
const_shts.ProjConst(end,:)=[];
for k=1:height(const_shts.ProjConst)
    const(k,:)=horzcat(const_shts.ProjConst(k,:),const_shts.(string(const_shts.ProjConst.tracker(k))));
end
for k=1:height(const)
    if isempty(const.proj_folder{k})
        const.proj_folder{k}=outpath;
    end
    if k==1
        folderpath = string(const.proj_folder{k});
    else
        folderpath = string(const.proj_folder{k});    filelist = dir(folderpath);
        if strcmpi(folderpath,'last')==1
            folderpath = string(const.proj_folder{k-1});
            const.proj_folder{k}=const.proj_folder{k-1};
        else
            folderpath = string(const.proj_folder{k});
        end
    end
    if k==1
        fptopo = append(folderpath,"/Topo/");
        filelist   = dir(fptopo);
        name       = {filelist.name};
        pat        =[".csv", ".xlsx"];
        name       = name(endsWith(name, pat, 'IgnoreCase',true));   % csv or xlsx files
        const.topo{k}=append(fptopo,char(name));
    else
        if strcmpi(string(const.topo(k)),"last")==1
            const.topo(k)=const.topo(k-1);
        else
            fptopo = append(folderpath,"/Topo/"); 
            filelist   = dir(fptopo);
            name       = {filelist.name};
            pat        =[".csv", ".xlsx"];
            name       = name(endsWith(name, pat, 'IgnoreCase',true));   % csv or xlsx files
            const.topo{k}=append(fptopo,char(name));
        end
    end
    filelist = dir(folderpath);
    name     = {filelist.name};
    npsname  = name(startsWith(name, 'NPS', 'IgnoreCase', true));
    spnname  = name(startsWith(name, 'SPANS', 'IgnoreCase', true));
    fpNPS = append(folderpath,'/',char(npsname)); 
    fpSPAN = append(folderpath,'/',char(spnname)); 
    const.keypile{k}=fpNPS;
    const.spans{k}=fpSPAN;
    const.cott_bins{k}=readtable(baseFileName,'UseExcel',true,'Sheet',const.cott_bins{k});
    if contains(const.tracker{k},'ojjo','IgnoreCase',true)
        const.geo_inputs{k}=readtable(baseFileName,'UseExcel',true,'Sheet',const.geo_inputs{k});
        const.geo_inputs{k}(:,3)=[];
        const.geo_inputs{k}=rows2vars(const.geo_inputs{k},"VariableNamesSource","GeoInputs");
        const.geo_inputs{k}(:,1)=[];
    end
    const.leg_factory{k}=readtable(baseFileName,'UseExcel',true,'Sheet',const.pile_factory{k});
    const.pile_bins{k}=readtable(baseFileName,'UseExcel',true,'Sheet',const.pile_binning{k});
    const.email{k}=readtable(baseFileName,'UseExcel',true,'Sheet','email');
    const.email{k}=rmmissing(const.email{k});
    bincheck=const.cott_bins{k};
    if strcmpi(string(const.grading_bins(k)),'all')==1
        minwp=bincheck.bin_start(1);
        maxwp=bincheck.bin_finish(end);
        const.min_wp(k)=minwp;
        const.max_wp(k)=maxwp;
    else
        minwp=bincheck.bin_start(1);
        const.min_wp(k)=minwp;
        if iscellstr(const.grading_bins(k))==1
            maxwp=bincheck.bin_finish(str2num(const.grading_bins{k}));
            const.max_wp(k)=maxwp;
        else
            maxwp=bincheck.bin_finish(const.grading_bins(k));
            const.max_wp(k)=maxwp;
        end
    end
    const.truss_leg_angle(k)=20;
    const.wp_slope_range(k)=60;
    if strcmp(const.tracker{k},'NEV')==1
        const.Nevados(k)=1;
    else
        const.Nevados(k)=0;
    end
    %const.outpath{k}=outpath;
    %const.t{k}=string(datetime("today"));
    mkdir(append(outpath,'/out_', string(datetime("today")),'/figures'));
    mkdir(append(outpath,'/out_', string(datetime("today")),'/data'));
    mkdir(append(outpath,'/out_', string(datetime("today")),'/rollups'));
    mkdir(append(outpath,'/out_', string(datetime("today")),'/grading'));
    const.fpath{k}=append(outpath,'/out_', string(datetime("today")),'/figures');
    const.dpath{k}=append(outpath,'/out_', string(datetime("today")),'/data');
    const.rupath{k}=append(outpath,'/out_', string(datetime("today")),'/rollups');
    const.gpath{k}=append(outpath,'/out_', string(datetime("today")),'/grading');
            %% Constants Setting - Just putting these here so they're not in the file
        if contains(const.tracker{k},'ojjo','IgnoreCase',true)
            const.trusscalc(k)=1; const.perp(k)=1; const.thirdleg(k)=0; %these are all ojjo requirements
            const.wp_slope_inrange=const.min_wp; const.wp_slope_outrange=const.min_wp; %this is an ojjo only requirement
        else
            const.trusscalc(k)=0; const.perp(k)=0; const.thirdleg(k)=0; %these are all ojjo requirements
            const.wp_slope_inrange=const.min_wp; const.wp_slope_outrange=const.min_wp; %this is an ojjo only requirement
        end
        const.writefiles(k)=1; const.simple_out(k)=1;
        %const.outpath{k}=[char(const.outpath{k}) '/' char(const.project{k}) '_out_' char(const.t{k})];% '/' char(const.project)];
end
    gsx(const.bin_x==2 & const.gsmooth==1)=3;
    gsx(const.bin_x==2 & const.gsmooth==2)=6;
    gsx(const.bin_x==2 & const.gsmooth==3)=12;
    gsx(const.bin_x==5 & const.gsmooth==1)=2;
    gsx(const.bin_x==5 & const.gsmooth==2)=4;
    gsx(const.bin_x==5 & const.gsmooth==3)=6;
    gsx(const.bin_x==10 & const.gsmooth==1)=1;
    gsx(const.bin_x==10 & const.gsmooth==2)=2;
    gsx(const.bin_x==10 & const.gsmooth==3)=3;
    gsy(const.bin_y==2 & const.gsmooth==1)=4;
    gsy(const.bin_y==2 & const.gsmooth==2)=8;
    gsy(const.bin_y==2 & const.gsmooth==3)=16;
    gsy(const.bin_y==5 & const.gsmooth==1)=2;
    gsy(const.bin_y==5 & const.gsmooth==2)=5;
    gsy(const.bin_y==5 & const.gsmooth==3)=7;
    gsy(const.bin_y==10 & const.gsmooth==1)=1;
    gsy(const.bin_y==10 & const.gsmooth==2)=2;
    gsy(const.bin_y==10 & const.gsmooth==3)=3;
    const.gsx=gsx';
    const.gsy=gsy';
const_all=const;
end