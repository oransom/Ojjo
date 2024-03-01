tic
warning ('off','all')
%uncomment following two lines for ojjo
% [FileName,PathName] = uigetfile('*.xlsx','Select the configuration file','MultiSelect','off'); %ojjo
% fp_cat=strcat(PathName,FileName); %cat path and file ojjo
% [const_all] = a_2_constants_config_6(fp_cat); % ojjo
%config_file='/Users/osop/Documents/ORCaS/Osc6_testproj/AUI_MtPoso_rev9SitePlan_022224/Constants_V6.0.xlsx';
%config_file='/Users/osop/Documents/ORCaS/Osc6_testproj/Dual String 1/Constants_V6.0.xlsx';
config_file='/Users/osop/Documents/ORCaS/Osc6_testproj/Zone1/BF2.xlsx';
%config_file='/Users/osop/Documents/ORCaS/Osc6_testproj/TheRyanCo_AVI-Fairfield_022224_V2/AVI_6.0.xlsx';
[const_all] = a_2_constants_config_6(config_file); %change to fp_cat for ojjo
for s=1:height(const_all.project)
    skips = [];
    if ~ismember(s,skips)
        fprintf('*************************************************** \n');
        fprintf('*******************~~~OSCaR6.0~~~****************** \n');
        fprintf('%s \n',string(const_all.project(s)));
        clearvars -except config_file s const_all surface
        if s>1
            if ~matches(string(const_all.topo(s)),string(const_all.topo(s-1)))
                clear surface
            end
        end
        close all
        const=const_all(s,:); %make sure you check hard coded consts below
        if ~contains(const.tracker,'Ojjo') %ojjo lock
            fprintf('***********************$$$$$$********************** \n');
            fprintf('********************PAY UP SUCKA******************* \n');
            return
        end
        %% create surface
        if exist('surface','var')==0
            [surface]=a_3_surfcreate_6(const);
        end
        fprintf('Creating topo complete %3.0f s \n',toc);

        %% northern pile and span excel files
        [drow, dpile, const] = a_4_npspans_6(const);

        if const.surftrim==1
            [surface]=a_5_surfacetrimmer_6(dpile, surface);
        end
        
        if const.np_trim==1
            [dpile] = a_6_pointstrimmer_6(dpile,surface);
        end

        %% Uncomment to check for pile/surface mismatch
        %a_7_surfcheck_6(const, drow, dpile, surface)
        
        %% Calculate all row fits
        [drow] = a_8_rowfits_6(const, surface, drow); %row best fit(s)

        %% pile location calculations: Fit, slope limit, flood, flipex, nsfit, motor align
        switch char(const.tracker)
            case {'ATI','NXT','Ojjo_ATI','Ojjo_NXT'} %straight torque tube solution
                [drow,dpile]=ab_1_plc_main_6(const,surface,dpile,drow);
            case {'XTR', 'XTR1p5'}
                [drow,dpile]=ac_1_xtrplc_main_6(const,surface,dpile,drow);
            case {'NEV'}
                [drow,dpile]=ad_1_nevplc_1_main(const,surface,dpile,drow);
            case {'FTC'}
                [drow,dpile]=xx_1_ftcplc_main_6(const,surface,dpile,drow);
        end
        fprintf('Pile location calcs complete %3.2f s \n', toc);
        fprintf('starting grading %3.2f s \n', toc);
        if strcmpi(const.gmethod,'none')
            [surface] = ae_1_grading_6(const,surface,dpile);
        else
            [surface, dpile] = ae_1_grading_6p1(const,surface,dpile,drow);
        end
        [surface] = ae_2_grading_extents_6(surface);
        if const.solve_postgrade==1
            [drow, dpile] = af_1_postgrading_6p1(const, surface, drow, dpile);
        end
        fprintf('finished grading %3.2f s \n', toc);

        fprintf('starting plotting %3.2f s \n', toc);
        [plots] = ah_1_planar_plot_6(const, surface, drow, dpile);
        fprintf('planar plot done %3.2f s \n', toc);
        [plots] = ah_2_slope_plot_6(const, surface, drow, dpile, plots);
        fprintf('slope plot done %3.2f s \n', toc);
        [plots] = ah_3_grading_plot_6(const, surface, plots);
        fprintf('grading plot done %3.2f s \n', toc);
        [plots] = ah_4_edge_plot_6(const, surface, drow, dpile, plots);
        fprintf('edge plot done %3.2f s \n', toc);

        %emailer2(const, kpv, pdfout, mname, pd_out, rp_loc, pbt) %changed to fname from pdfout for apple osx testing

        %% clean it up
        xx_cleanup_6;

        fprintf('Full program run complete %3.2f s \n', toc);
        fprintf('*************************************************** \n');

    end
end

