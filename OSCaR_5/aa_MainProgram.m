tic
warning ('off','all')
%uncomment following two lines for ojjo
%[FileName,PathName] = uigetfile('*.xlsx','Select the configuration file','MultiSelect','off');
%fp_cat=strcat(PathName,FileName); %cat path and file
config_file='/Volumes/ORCaS/FTC/River Trail/Pioneer/Constants_V5.04_RT_Pio.xlsx';
[const_all] = constants_config(config_file); %change to fp_cat for ojjo
for s=1:height(const_all.project)
    skips = [];
    if ~ismember(s,skips)
        fprintf('************ %s ************\n',string(const_all.project(s)));
        clearvars -except config_file s const_all perm_const Flood_s surface F_og%Flood
        if s>1
            if matches(string(const_all.topo(s)),string(const_all.topo(s-1))) == 0
                clear surface
                clear F_og
                clear Flood_s
            end
        end
        close all
        const=const_all(s,:); %make sure you check hard coded consts below
        const.simple_out=1; %lock in simple out
        const.writefiles=1; %lock in writing files
        %% Constants Setting - Just putting these here so they're not in the file
        ojjo_lock=0;
        if contains(const.tracker,'ojjo','IgnoreCase',true)
            const.trusscalc=1; const.perp=1; const.thirdleg=0; %these are all ojjo requirements
            const.wp_slope_inrange=const.min_wp; const.wp_slope_outrange=const.min_wp; %this is an ojjo only requirement
        else
            const.trusscalc=0; const.perp=0; const.thirdleg=0; %these are all ojjo requirements
            const.wp_slope_inrange=const.min_wp; const.wp_slope_outrange=const.min_wp; %this is an ojjo only requirement
        end
        const.sect_name=const.project; %legacy - should take out
        const.outpath=[char(const.outpath{1}) '/' char(const.project{1}) '_out_' char(const.t{1})];% '/' char(const.sect_name)];
        %%
        if exist('surface','var')==0
            gtopo=readtable(string(const.topo));
            gtopo=rmmissing(gtopo);
            if strcmpi(const.flood,'na')==0
                if isnumeric(const.flood)==0
                    flood=readtable(string(const.flood)); %flood
                end
            end
        end
        %% create surface
        if exist('surface','var')==0
            surface.gtopox=gtopo.Var1;
            surface.gtopoy=gtopo.Var2;
            surface.gtopoz=gtopo.Var3;
            F_og=scatteredInterpolant(surface.gtopox,surface.gtopoy,surface.gtopoz,'natural','none');
            clear surface;
            %% flood stuff
            if strcmpi(const.flood,'na')==0
                if isnumeric(const.flood)==0
                    Flood_s=scatteredInterpolant(flood.Var1,flood.Var2,flood.Var3,'natural','none'); %flood
                else
                    Flood_s=F_og;
                end
            else
                Flood_s=0;
            end
            [surface.xq,surface.yq]=meshgrid(min(gtopo.Var1):const.bin_x:max(gtopo.Var1), min(gtopo.Var2):const.bin_y:max(gtopo.Var2));
            clear gtopo;
            surface.zog=F_og(surface.xq,surface.yq);
            % surface.flood=Flood_s(surface.xq,surface.yq);
            % surface.flood(surface.flood<surface.zog)=surface.zog(surface.flood<surface.zog);
            % Flood_s=scatteredInterpolant(reshape(surface.xq,[],1),reshape(surface.yq,[],1),reshape(surface.flood,[],1),'natural','none');
            % surface.flood=Flood_s(surface.xq,surface.yq);
            % foo=1;
            %% flood stuff
            if ~strcmpi(const.flood,'na')
                if ~isnumeric(const.flood) & contains(const.flood,'.csv')
                    surface.flood=Flood_s(surface.xq,surface.yq); %flood
                elseif isnumeric(const.flood) %if a numeric value is given set everything below that value to flood depth given
                    sfoo=surface.zog<const.flood;
                    sfoo=sfoo.*surface.zog;
                    sfoo(sfoo==0)=NaN;
                    surface.flood=const.flood-sfoo;
                    Flood_s.Values=((Flood_s.Values<const.flood).*const.flood)-F_og.Values.*(Flood_s.Values<const.flood);
                end
                surface.floodmask=isnan(surface.flood);
                if mean(mean(surface.flood,'omitnan'),'omitnan')<mean(mean(surface.zog,'omitnan'),'omitnan')
                    surface.flood(isnan(surface.flood))=0;
                    surface.flood(surface.flood<0)=0;
                    surface.flood=surface.flood+surface.zog;
                    Flood_s=griddedInterpolant(surface.xq',surface.yq',surface.flood','nearest','none'); %flood
                end
            end
            clear surface.gtopox
            clear surface.gtopoy
            clear surface.gtopoz
        end

        if s==1
            create_topo=toc;
            fprintf('Creating topo complete %3.0f s \n',create_topo);
        end
        %% load in all northern pile excel files
        Folder = string(const.keypile);
        if isfolder(Folder) ~= 1
            Message = sprintf('Error: The following folder does not exist:\n%s', Folder);
            uiwait(warndlg(Message));
            return;
        end
        filePattern = fullfile(Folder, '*.xlsx');
        pileFiles   = dir(filePattern);
        for k = 1:length(pileFiles)
            baseFileName = pileFiles(k).name;
            fullFileName = fullfile(Folder, baseFileName);
            kpv{k}=erase(pileFiles(k).name, ".xlsx");
            try
                dat.(kpv{k}) = readtable(fullFileName);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:spreadsheet:book:fileOpen'))
                    msg = 'Close all dependent project files';
                    causeException = MException('MATLAB:spreadsheet:book:fileOpen',msg);
                    ME = addCause(ME,causeException);
                end
                rethrow(ME)
            end
        end
        if sum(contains(fieldnames(dat.(kpv{1})),'row_length'))==0
            %% load in all span types
            Folder = string(const.spans);
            if isfolder(Folder) ~= 1
                Message = sprintf('Error: The following folder does not exist:\n%s', Folder);
                uiwait(warndlg(Message));
                return;
            end
            filePattern = fullfile(Folder, '*.xlsx');
            spanFiles   = dir(filePattern);
            for k = 1:length(spanFiles)
                baseFileName = spanFiles(k).name;
                fullFileName = fullfile(Folder, baseFileName);
                sv{k}=erase(spanFiles(k).name, ".xlsx");
                try
                    span.(sv{k}) = readtable(fullFileName);
                catch ME
                    if (strcmp(ME.identifier,'MATLAB:spreadsheet:book:fileOpen'))
                        msg = 'Close all dependent project files';
                        causeException = MException('MATLAB:spreadsheet:book:fileOpen',msg);
                        ME = addCause(ME,causeException);
                    end
                    rethrow(ME)
                end
            end
        else
            for i=1:length(kpv)
                [foo,fooi]=unique(dat.(kpv{i}).span);
                for j=1:length(foo)
                    span.(foo{j}).s=foo{j};
                    if span.(foo{j}).s<140
                        span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/6:dat.(kpv{i}).row_length(fooi(j))];
                        span.(foo{j}).kps=span.(foo{j}).kps';
                        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                        for k=1:length(span.(foo{j}).kps)
                            span.(foo{j}).tl{k,1}=Alphabet(k);
                            if k~=4
                                span.(foo{j}).tt{k,1}='std';
                            else
                                span.(foo{j}).tt{k,1}='mtr';
                            end
                        end
                    elseif span.(foo{j}).s>140 && span.(foo{j}).s<250
                        span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/10:dat.(kpv{i}).row_length(fooi(j))];
                        span.(foo{j}).kps=span.(foo{j}).kps';
                        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                        for k=1:length(span.(foo{j}).kps)
                            span.(foo{j}).tl{k,1}=Alphabet(k);
                            if k~=6
                                span.(foo{j}).tt{k,1}='std';
                            else
                                span.(foo{j}).tt{k,1}='mtr';
                            end
                        end
                    elseif span.(foo{j}).s>250 && span.(foo{j}).s<350
                        span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/12:dat.(kpv{i}).row_length(fooi(j))];
                        span.(foo{j}).kps=span.(foo{j}).kps';
                        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                        for k=1:length(span.(foo{j}).kps)
                            span.(foo{j}).tl{k,1}=Alphabet(k);
                            if k~=7
                                span.(foo{j}).tt{k,1}='std';
                            else
                                span.(foo{j}).tt{k,1}='mtr';
                            end
                        end
                    else
                        span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/14:dat.(kpv{i}).row_length(fooi(j))];
                        span.(foo{j}).kps=span.(foo{j}).kps';
                        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                        for k=1:length(span.(foo{j}).kps)
                            span.(foo{j}).tl{k,1}=Alphabet(k);
                            if k~=8
                                span.(foo{j}).tt{k,1}='std';
                            else
                                span.(foo{j}).tt{k,1}='mtr';
                            end
                        end
                    end
                end
            end
        end


        % call to the pile place function, sends F_og (surface), kpv (list of input
        % blocks), dat (raw north pile data), span (raw span data), constant data
        % returns calculated pile location (cpl) and some pile meta data (pmd), truss
        % data (trs) and out of spec data (oos_data)
        % % %
        % figure
        % %surf(surface.xq,surface.yq,surface.zog)%_trimmed)
        % surf(surface.xq,surface.yq,surface.flood);
        % hold on
        % contour3(surface.xq,surface.yq,surface.zog,500)
        % shading interp
        % alpha(0.5)
        % view(2)
        % hold on
        % for i=1:numel(kpv)
        % scatter(dat.(kpv{i}).tpx,dat.(kpv{i}).tpy,15,'filled')
        % scatter(dat.(kpv{i}).tpx,dat.(kpv{i}).tpy-dat.(kpv{i}).row_length,15,'filled')
        % end
        % daspect([1 1 1])


        %topo->NP trimmer - eliminate any NP that aren't within topo
        %footprint
        const.surftrim=1;
        surface.untrimmed=surface.zog;
        if const.surftrim==1
            [surface]=surfacetrimmer(kpv,dat,span,const,surface);
        end
        if const.np_trim==1
            [dat] = NPtrimmer(kpv,dat,span,const,surface);
        end


        try
            if const.Nevados==1
                [cpl, uc, trs, pmd, oos_data, ns_ugly] = pilelocationcalc_Nev(F_og,Flood_s,kpv,dat,span,const,surface); %Flood, Nevados
            elseif strcmp(const.tracker,'XTR')==1
                [cpl, uc, trs, pmd, oos_data, ns_ugly] = pilelocationcalc_XTR_1p1(F_og,Flood_s,kpv,dat,span,const,surface);
            elseif strcmp(const.tracker,'XTR1p5')==1
                [cpl, uc, trs, pmd, oos_data, ns_ugly] = pilelocationcalc_XTR1p5(F_og,Flood_s,kpv,dat,span,const,surface);
            else
                [cpl, uc, trs, pmd, oos_data, ns_ugly] = pilelocationcalc(F_og,Flood_s,kpv,dat,span,const,surface); %Flood
            end
        catch ME
            if (strcmp(ME.identifier,'MATLAB:nonExistentField'))
                msg = 'Most likely you are calling a span table from keypile file that does not exist';
                causeException = MException('MATLAB:spreadsheet:book:fileOpen',msg);
                ME = addCause(ME,causeException);
            end
            rethrow(ME)
        end
        plc=toc;
        fprintf('Pile location calcs complete %3.2f s \n', plc);

        [po, to, pmd, lengthst] = pile_truss_output(kpv, const, cpl, trs, pmd);
        for i=1:length(kpv)
            if i==1
                oos_master=[po.(kpv{i}).os_pl];
            else
                oos_master=[oos_master;po.(kpv{i}).os_pl];
            end
        end
        pto=toc;
        fprintf('Pile truss output complete %3.2f s \n', pto);
        % create pinning table
        [pin_out] = pin_table(kpv, pmd, po);  %CAN WE MOVE THIS DOWN?
        % create construction docs table
        [con_out] = construction_docs(kpv, pmd, po); %CAN WE MOVE THIS DOWN?
        % create pier plot plan
        [ppp_out] = pierplot(kpv, pmd, po, const); %CAN WE MOVE THIS DOWN?
        ppo=toc;
        fprintf('Pier plot complete %3.2f s \n', ppo);
        try
            if const.Nevados==1
                [pd_out, ot, pileTable] = piledata_Nev(kpv, pmd, po, to, oos_data, const);
            else
                [pd_out, lfout, ot, pileTable] = piledata(kpv, pmd, po, to, oos_data, const);
            end
        catch ME
            if (strcmp(ME.identifier,'MATLAB:struct2table:UnequalFieldLengths'))
                msg = 'Most likely you have a truss type naming issue in your span tables not matching one of Standard, Standard Motor, Heavy, Heavy Motor, or you have duplicate sections in a NP file';
                causeException = MException('MATLAB:spreadsheet:book:fileOpen',msg);
                ME = addCause(ME,causeException);
            end
            rethrow(ME)
        end
        %% search for nearest neighbors out of spec - call function oos
        if matches(string(const.tracker),"NXT") || contains(string(const.tracker),"XTR",'IgnoreCase',true)==1
            [rnbr_table, mover_table, fliped] = rowmover_test4_nxt_edg(const, ot);
            [rnbr_table, mover_table, flipex] = rowmover_test4_nxt_ext(const, ot);
        else
            [rnbr_table, mover_table, flipex] = rowmover_test4(const, ot);
            fliped=[];
        end
        % this loop plugs what has been moved (mover_table.raiserow) back into the con_out and ppp construction documents
        % it overwrites values with E/W raiserow values so only one set of values exists
        for i=1:length(kpv)
            for j=1:length(con_out.(kpv{i}).bsr)
                if find(strcmp(con_out.(kpv{i}).bsr(j),mover_table.bsr))>0
                    foo = find(strcmp(con_out.(kpv{i}).bsr(j),mover_table.bsr));
                    con_out.(kpv{i}).WP_cm_final(j)=con_out.(kpv{i}).WP_cm(j)+mover_table.raiserow(foo)*30.48; %raise row by cm
                    ppp_out.(kpv{i}).reveal_ht_ft(j)=ppp_out.(kpv{i}).reveal_ht_ft(j)+mover_table.raiserow(foo);
                    ppp_out.(kpv{i}).reveal_ht_cm(j)=(ppp_out.(kpv{i}).reveal_ht_ft(j))*30.48;
                    ppp_out.(kpv{i}).top_pile_elev_ft(j)=ppp_out.(kpv{i}).top_pile_elev_ft(j)+mover_table.raiserow(foo);
                else
                    con_out.(kpv{i}).WP_cm_final(j)=con_out.(kpv{i}).WP_cm(j);
                    ppp_out.(kpv{i}).reveal_ht_ft(j)=ppp_out.(kpv{i}).reveal_ht_ft(j);
                    ppp_out.(kpv{i}).reveal_ht_cm(j)=ppp_out.(kpv{i}).reveal_ht_cm(j);
                    ppp_out.(kpv{i}).top_pile_elev_ft(j)=ppp_out.(kpv{i}).top_pile_elev_ft(j);
                end
            end
            con_out.(kpv{i}).TD3=con_out.(kpv{i}).WP_cm_final>=245;
            td3bsr=unique(con_out.(kpv{i}).bsr(con_out.(kpv{i}).TD3==1),'stable');
            if ~isempty(td3bsr)
                for j=1:length(td3bsr)
                    mask=strcmp(con_out.(kpv{i}).bsr,td3bsr(j));
                    con_out.(kpv{i}).TD3(mask)=1;
                end
            end
        end
        %next line needs to be commented for non-ojjo uses
        con_out.(kpv{i})=[con_out.(kpv{i})(:,1:9) con_out.(kpv{i})(:,11) con_out.(kpv{i})(:,10)];
        cco=toc;
        fprintf('Construction docs complete %3.2f s \n', cco);
        %pointfile
        [top_out] = pile_top(kpv, pmd, po, ppp_out);
        %output files loop
        if contains(const.tracker,'ojjo','IgnoreCase',true) && const.writefiles == 1
            %mkdir('output');
            for i=1:length(kpv)
                mkdir([const.outpath '/ojjo/' kpv{i} '_output']);
                dlmwrite([const.outpath '/ojjo/' kpv{i} '_output/' 'outspec_pile.csv'],po.(kpv{i}).os_pl, 'delimiter', ',', 'precision', '%0.4f');
                dlmwrite([const.outpath '/ojjo/' kpv{i} '_output/' 'outspec_nslope.csv'],po.(kpv{i}).os_ns, 'delimiter', ',', 'precision', '%0.4f');
                dlmwrite([const.outpath '/ojjo/' kpv{i} '_output/' 'outspec_pslope.csv'],po.(kpv{i}).os_ps, 'delimiter', ',', 'precision', '%0.4f');
            end
        end

        if const.writefiles == 1
            writetable(mover_table,[const.outpath '/flipex/' char(const.sect_name) ' move-rows.csv'])
        end
        if const.writefiles == 1
            mkdir([const.outpath '/output']);
            if contains(string(const.project),"genpile","IgnoreCase",true)==1
                mkdir([const.outpath '/output/.xtemp/']);
                dlmwrite([const.outpath '/output/.xtemp/' kpv{i} '_OoS_Pile.csv'],po.(kpv{i}).os_pl, 'delimiter', ',', 'precision', '%0.4f');
            end
            for i=1:length(kpv)
                dlmwrite([const.outpath '/output/' kpv{i} '_OoS_Pile.csv'],po.(kpv{i}).os_pl, 'delimiter', ',', 'precision', '%0.4f');
            end
        end

        %grading function and plotter
        grad=toc;
        fprintf('starting grading %3.2f s \n', grad);
        [total_grad,grd_ext,grd_area,graded_xyz,gname,gnamepdf] = bfgrade3(const, kpv, po, surface, ppp_out, F_og); %switching to bfgrade2 for testing 05/22/23
        if grd_area<1 || const.solve_postgrade==0
            xyzf_int=F_og;
            if const.writefiles==1
                mkdir([const.outpath '/output/']);
                for i=1:length(kpv)
                    writetable(pd_out.(kpv{i}),[const.outpath '/output/' 'Pile_Data_all.csv']);
                end
            end
        else
            if contains(const.tracker,'ojjo','IgnoreCase',true)
                [pd_out,con_out,ppp_out,xyzf_int,top_out,pin_out] = legdata_postgrade(kpv, graded_xyz, pd_out, const, con_out, ppp_out, top_out, pin_out); %graded_xyz is required
            else
                [pd_out,piletable_upr,ppp_out,xyzf_int] = postdata_postgrade(kpv, graded_xyz, pd_out, const, ppp_out); %graded_xyz is required
            end
        end
        graddone=toc;
        fprintf('finished grading %3.2f s \n', graddone);
        %plotpile
        plotpile(const, surface, kpv, po, cpl, flipex, ns_ugly, pd_out)
        plotpile_flood(const, surface, kpv, pd_out)
        pltpile=toc;
        fprintf('Plotpile complete %3.2f s \n', pltpile);

        if const.Nevados ~=1
            [master_dat]=mastercombiner_pile(kpv, ot, const, pd_out, lfout, pileTable, mover_table, con_out, ppp_out, total_grad, flipex, ns_ugly, oos_master,grd_area);
        else
            lfout=[];
            [master_dat]=mastercombiner_pile(kpv, ot, const, pd_out, lfout, pileTable, mover_table, con_out, ppp_out, total_grad, flipex, ns_ugly, oos_master,grd_area);
        end
        if const.Nevados ~= 1 && ~contains(const.tracker,'ojjo','IgnoreCase',true)
            [rpcs_dat, rp_loc]=rpcs_table(kpv, ot, const, pd_out, pileTable, mover_table, con_out, ppp_out, total_grad, flipex, ns_ugly, oos_master);
        else
            rp_loc=[];
        end
        if const.Nevados==1
            nev_plot(kpv,pd_out,F_og,surface,ppp_out, const)
        end

        plot=toc;
        fprintf('starting plotting %3.2f s \n', plot);
        [fname, fnamepdf, pbt] = rpcs_plot3(kpv,pd_out,F_og,xyzf_int,surface,ppp_out, const, grd_ext, fliped, flipex);
        [sname, snamepdf] = rpcs_plot_slope_new(kpv,pd_out,F_og,surface,ppp_out, const,grd_ext, fliped, flipex);
        %rpcs_plot_buildable(kpv,pd_out,F_og,surface,ppp_out, const, grd_ext, flipex)
        plotdone=toc;
        fprintf('finished plotting %3.2f s \n', plotdone);
        %rpcs_fhm(kpv,pd_out,F_og,surface,ppp_out, const, grd_ext)
        if ~contains(const.tracker,'ojjo','IgnoreCase',true)
            [mname] = sum_of_sums(const, kpv, pd_out);
        else
            mname=[];
        end
        full_loop=toc;
        fprintf('Full loop complete %3.2f s \n', full_loop);
        if const.Nevados~=1 && contains(string(const.project),"genpile","IgnoreCase",true)~=1
            rmdir([const.outpath '/output/.xtemp/'],'s');
        end
        %gname sname fname
        if ispc
            if const.Nevados~=1
                pdfout=strrep(fnamepdf,'_Planar_','_MERGED_PDF_');
                pdfout=strrep(pdfout,'pptx','pdf');
            else
                pdfout=strrep(fnamepdf,'_Eng_Study_','_MERGED_PDF_');
                pdfout=strrep(pdfout,'pptx','pdf');
            end
            if ~isempty(gnamepdf)
                mergePdfs([cellstr(fnamepdf),cellstr(gnamepdf),cellstr(snamepdf)],pdfout)
            else
                mergePdfs([cellstr(fnamepdf),cellstr(snamepdf)],pdfout)
            end
            delete(snamepdf);
            delete(gnamepdf);
            delete(fnamepdf);
        else
            % pptxout = strrep(fname,'_Planar_','_MERGED_PPT_');
            % mpptx    = exportToPPTX('', ...
            %     'Dimensions',[17 11], ...
            %     'Title','Merged OSCaR Output', ...
            %     'Author','ORCaS Inc., 2024', ...
            %     'Subject','Automatically generated documentation', ...
            %     'Comments','This file has been automatically generated by ORCaS_OSCaR_2024');
            % gpptx    = exportToPPTX(gname);
            % spptx    = exportToPPTX(sname);
            % fpptx    = exportToPPTX(fname);
            % %mslide1     = mpptx.addSlide('Master',fpptx);
            % %mslide1  = gpptx;
            % newFile = mpptx.save(pptxout);
            % %newFile = fpptx.save(pptxout);
            % foo=1;
        end
        %emailer2(const, kpv, pdfout, mname, pd_out, rp_loc, pbt) %changed to fname from pdfout for apple osx testing
    end
end
full_program=toc;
fprintf('Full program run complete %3.2f s \n', full_program);
fprintf('*************************************************** \n');