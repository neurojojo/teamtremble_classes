classdef displacedMouse < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    % To see how the AUC is calculated check line 356 (auc_rect) variable
    % which produces the calibration factor for the auc score
    
    properties
        mouse
        displacement
        displacement_csplus
        displacement_csminus
        csplus_index
        csminus_index
        ingress_index_csplus
        ingress_index_csminus
        displacement_csplus_filtered
        displacement_csminus_filtered
        fitlines_CSplus
        fitlines_CSminus
        corrected_displacement_csplus
        corrected_displacement_csminus
        h5_datasets
        cwt_csplus
        cwt_csminus
        ingress_stats_CSplus
        ingress_stats_CSminus
        aucs_fraction
        trembleScore_csplus
        trembleScore_csminus
        ingressScore_csplus
        ingressScore_csminus
        datadir
        stimulus_order
        keyTable
        keyString
        trialStruct
        trialTable
        ingressTable
        ingressBuffer = 1000;
        cwt_window = [1:150];
        
        t0 = 60000;
        tend = 200000;
        dsFactor = 10;
        
    end
    
    properties( Access = private )
        
        displacement_filename = 'DisplacementInOdor.mat';
        stimulus_filename = 'stimuli.mat';
        % FOR CWT SPACE ISSUE %
        % constrain the size to this value %
        w_limits_x = [200:400];
        
        threshold_min = 0.005; % Used in getTrembleScore
        
        fs = 1000; % Downsampled sampling frequency
        trials = 20;
        time_vector = (0:14000)/1000;
        
        Nframes = [];
        
        d = 400; % An internal parameter that is calibrated to the duration before and after an ingress
        % Use this in the fitLine subfunction

        
    end
    
    methods
        
        function obj = displacedMouse(mouseid, matfile_location)
            
            obj.mouse = mouseid;
            obj.datadir = matfile_location;
            obj.Nframes = numel( obj.time_vector );
            
        end
        
        function loadDisplacements( obj )
            fprintf('Loading displacement for %s\n', obj.mouse );
            
            displacementfile = fullfile( obj.datadir, obj.displacement_filename );
            
            displacement_Mat = load( displacementfile );
            obj.displacement = displacement_Mat.displacement(:,[obj.t0:obj.dsFactor:obj.tend],:);
            
            stimulusfile = fullfile( obj.datadir, obj.stimulus_filename );
            
            stimulus_Mat = load( stimulusfile );
            obj.stimulus_order = stimulus_Mat.stimuli;
            
            if ~isempty(obj.csplus_index)
                obj.displacement_csplus = reshape( permute( obj.displacement(:,:, ...
                                obj.csplus_index(isnan(obj.csplus_index)==0)  ), [2,1,3] ), size(obj.displacement,2), ...
                                prod( [size(obj.displacement,1), sum(isnan(obj.csplus_index)==0)] ) )';
                            % The index of the matrix obj.displacement will
                            % specify the CS type %
                fprintf('Filtering\n');
                %obj.displacement_csplus_filtered = cell2mat( ...
                %    arrayfun(@(x) filtfilt( obj.butter_b, obj.butter_a, obj.displacement_csplus(x,:) ), [1:obj.trials],...
                %    'UniformOutput', false, 'ErrorHandler', @(x,y) zeros(1,size(obj.displacement_csplus,2)) )' );
            end
            
            if ~isempty(obj.csminus_index)
                obj.displacement_csminus = reshape( permute( obj.displacement(:,:, ...
                                obj.csminus_index(isnan(obj.csminus_index)==0)  ), [2,1,3] ), size(obj.displacement,2), ...
                                prod( [size(obj.displacement,1), sum(isnan(obj.csminus_index)==0)] ) )';
                %obj.displacement_csminus_filtered = cell2mat( ...
                %    arrayfun(@(x) filtfilt( obj.butter_b, obj.butter_a, obj.displacement_csminus(x,:) ), [1:obj.trials],...
                %    'UniformOutput', false, 'ErrorHandler', @(x,y) zeros(1,size(obj.displacement_csminus,2)) )' );
            end
            
            if and(~isempty(obj.displacement_csminus),~isempty(obj.displacement_csplus))
                obj.displacement = [];
            end
            
        end
        
        function selectIngress(obj)

            csplus_locs = arrayfun(@(x) obj.changepts( obj.displacement_csplus(x,:), 'minpeakheight', 2.5 ), [1:size( obj.displacement_csplus ,1)], 'UniformOutput', false );
            csminus_locs = arrayfun(@(x) obj.changepts( obj.displacement_csminus(x,:), 'minpeakheight', 2.5 ), [1:size( obj.displacement_csminus ,1)], 'UniformOutput', false );

            % Obtaining the change points by manual inspection
            % (CS plus)

            keypressFxn1 = @(x) fprintf('%i', str2double(x) );

            f = figure('WindowKeyPressFcn', @(handle,event) keypressFxn1(event.Key), 'Tag', 'Fig1', 'WindowState','maximized' );

            ingressTable = table ( [1:size(csplus_locs,2)+size(csminus_locs,2)]',...
                [csplus_locs';csminus_locs'], [obj.displacement_csplus;obj.displacement_csminus],...
                [repmat({'csplus'},size(csplus_locs,2),1);repmat({'csminus'},size(csplus_locs,2),1)],...
                nan( size(csplus_locs,2)+size(csminus_locs,2), 1 ),...
                'VariableNames', {'Index','Locs','Displacement','Type','Changepoint'});
            idx = [1:size(obj.displacement_csplus,2)];
            
            for i=1:size(ingressTable,1)

                plot( ingressTable(i,:).Displacement ); ylim( [ min(ingressTable(i,:).Displacement), max(ingressTable(i,:).Displacement) ] );
                lines_ = ingressTable(i,:).Locs{1};
                title( sprintf('%s Trial %i',obj.mouse,i) );
                rowfun( @(x,y) text(y, 0.9*max(ingressTable(i,:).Displacement) ,sprintf('%i',x) ), table( [1:numel(lines_)]', lines_' ) );
                arrayfun( @(x) line([x,x],[-20,30],'color','k'), lines_ );
                
                tmp = str2num(cell2mat(inputdlg('Enter ingress line')));
                if ne(tmp,0); ingressTable(i,:).Changepoint = ingressTable(i,:).Locs{1}( tmp ); else; ingressTable(i,:).Changepoint = 10000; end
                
            end
            
            obj.ingressTable = ingressTable;
            fprintf('\nSaved ingress table for %s\n',obj.mouse);
            
        end
        
        % Set properties of this mouse %
        function associate_csplus_index(obj, array_ingress_times)
            obj.ingress_index_csplus = array_ingress_times;
        end
        
        function associate_csminus_index(obj, array_ingress_times)
            obj.ingress_index_csminus = array_ingress_times;
        end


        % Line-fitting %
        function myfitline = fitLine( obj, trial, trialType, varargin ) % Fits a piecewise linear function
            d = 1000;
            Nframes = obj.Nframes;
            
            fprintf('Computing for trial %i\n', trial );
            
            %%%%%%%%%%%%%%%%%%
            % Error handling %
            %%%%%%%%%%%%%%%%%%
            
            if and( strcmp(varargin{1},'coeffs'), size(obj.displacement_csplus,1) < trial)
                myfitline = nan(1,5);%'No data for this trial';
               return
            end

            if and( strcmp(varargin{1},'fitlines'), size(obj.displacement_csplus,1) < trial)
                myfitline = nan( 1, 14001 );
               return
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check for whether the trialType is CSplus or CSminus %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(trialType,'CSplus')
                Y = obj.displacement_csplus(trial,:);
                ing = obj.ingress_index_csplus(trial);
                % WE SET NO INGRESSES TO 10,000
                if eq(ing,10000); myfitline = smooth(Y,Nframes/10)'; return; end; % Create zero-vector for displacements lacking ingress
            else if strcmp(trialType,'CSminus')
                Y = obj.displacement_csminus(trial,:);
                ing = obj.ingress_index_csminus(trial);
                % WE SET NO INGRESSES TO 10,000
                if eq(ing,10000); myfitline = smooth(Y,Nframes/10)'; return; end; % Create zero-vector for displacements lacking ingress
                else
                    fprintf("Specify either 'CSplus' or 'CSminus' as a third argument\n");
                    myfitline = [];
                    return
                end
            end
               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find coefficients for the lines %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp( varargin{1}, 'coeffs' )

                Y_ = Y;
                Y_([1:ing-d]) = 0;
                Y_([ing-d+1:min( numel(Y), ing+d )]) = Y([ing-d+1:min( numel(Y), ing+d )]);
                Y_(ing+d+1:end) = mean(Y(ing+d+1:end));

                [pks_min,locs_min] = min( Y_ );
                [pks_max,locs_max] = max( Y_ );

                % Parabolic fit
                coeffs = polyfit( [0:locs_max-locs_min], Y_(locs_min:locs_max), 2 );
                    
                myfitline = struct();
                if locs_min > locs_max
                    fprintf('Error for %s on trial %i. Check ingress time!\n', obj.mouse, trial)
                    myfitline = nan(1,5);%'Low quality';
                else
                    myfitline = [locs_min,locs_max,coeffs];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create the fitted lines if coefficients exist        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp( varargin{1}, 'fitlines' )
                
                coeffs = obj.(sprintf('ingress_stats_%s',trialType)){trial};
                
                if any(isnan(coeffs)); myfitline = nan( 1, 14001 ); return; end
                
                [locs_min,locs_max] = deal( coeffs(1), coeffs(2) );
                coeffs = coeffs(3:end);
                
                Fit_Table = struct();
                Fit_Table.X_pre = [0:locs_min-1];
                Fit_Table.Y_pre = zeros( numel(1, Fit_Table.X_pre ), 1 );
                Fit_Table.X_ing = [locs_min:locs_max-1];
                Fit_Table.Y_ing = polyval( coeffs, [locs_min:locs_max-1]-locs_min )';
                Fit_Table.X_post = [locs_max+1:numel(Y)];
                Fit_Table.Y_post = mean(Y([locs_max+1:numel(Y)])) * ones( numel(Fit_Table.X_post), 1 );

                myfitline = [ Fit_Table.Y_pre; Fit_Table.Y_ing; Fit_Table.Y_post ]';
                
            end
            
            
        end
        
        function result = searchMouse( obj, query )
        
            output = 1+(regexp( obj.keyString, query ) - 1)/2;
            query_stripped = regexp(query,'[pm]','match');
            len_sequence = numel( query_stripped );
            
            % Retrieve based on csp or csm
            result = [];
            table_variables = {};
            if ~isempty(output)
                for i = output;
                    tmp_result = [];
                    for j = 0:len_sequence-1
                        tmp_result = [ tmp_result, obj.trialTable( i+j, : ).(query_stripped{j+1})];
                        table_variables{j+1} = sprintf('%s%i',query_stripped{j+1},j+1);
                    end
                    result = [result; tmp_result];
                end
            end
            result = array2table( result, 'VariableNames', table_variables );
        end
        
        function getTrembleScore( obj )
           
            fprintf('Computing tremble score for %s\n', obj.mouse );
            tic;
            No = 10;
            Nv = 48;
            d = 800;
            
            % CSminusPLUS
            fprintf('Computing cwt for %s\n',obj.mouse);
            
            if isempty( obj.cwt_csplus )
                w = arrayfun( @(trial) cwt( obj.corrected_displacement_csplus( trial , :), obj.fs, 'NumOctaves',No, 'VoicesPerOctave',Nv), [1:obj.trials], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN );
                close all;
                w = cellfun( @(this_w) this_w( obj.w_limits_x, : ), w, 'ErrorHandler', @(x,y) NaN, 'UniformOutput', false );
                w_conj = arrayfun( @(trial) w{trial}.*conj( w{trial} ), [1:obj.trials], 'UniformOutput', false );
                % REMOVE EXCESSIVE INFORMATION FROM CWT %
                w_conj = cellfun( @(this_cell) this_cell .* im2bw( this_cell,10^-5 ), w_conj, 'UniformOutput', false );
                obj.cwt_csplus = w_conj;
            end
            
            % Frequency band 200-350, time is 2000:8000 (2000 = cs on, 8000
            % = arbitrary end)
            ingressBuffer = 1000;
            
            obj.trembleScore_csplus.Mag = cell2mat( arrayfun( @(i) obj.max_of_grad_adj_img( obj.cwt_csplus{i}(:,[2000: obj.ingress_index_csplus(i)-ingressBuffer ])), [1:numel(obj.ingress_index_csplus)], 'ErrorHandler', @(x,y) NaN, 'UniformOutput', false ) );
            
            if numel(obj.trembleScore_csplus.Mag)<20;
                1
            end
            
            % CSminusMINUS
            if isempty( obj.cwt_csminus )
                w = arrayfun( @(trial) cwt( obj.corrected_displacement_csminus( trial , :), obj.fs, 'NumOctaves',No, 'VoicesPerOctave',Nv), [1:obj.trials], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN, 'UniformOutput', false );
                close all;
                w = cellfun( @(this_w) this_w( obj.w_limits_x, : ), w, 'ErrorHandler', @(x,y) NaN, 'UniformOutput', false );
                w_conj = arrayfun( @(trial) w{trial}.*conj( w{trial} ), [1:obj.trials], 'UniformOutput', false );
                % REMOVE EXCESSIVE INFORMATION FROM CWT %
                w_conj = cellfun( @(this_cell) this_cell .* im2bw( this_cell,10^-5 ), w_conj, 'UniformOutput', false );
                obj.cwt_csminus = w_conj;
            end
            
            % Frequency band 200-350, time is 2000:8000 (2000 = cs on, 8000
            % = arbitrary end)
            
            obj.trembleScore_csminus.Mag = cell2mat( arrayfun( @(i) obj.max_of_grad_adj_img( obj.cwt_csminus{i}(:,[2000: obj.ingress_index_csminus(i)-ingressBuffer ])), [1:numel(obj.ingress_index_csminus)], 'ErrorHandler', @(x,y) NaN, 'UniformOutput', false ) );
                        
            if numel(obj.trembleScore_csminus.Mag)<20;
                1
            end
            
            % Duration = Fraction of time between CS onset and the detected
            % ingress during the threshold tremble power of 0.001 is exceeded
            obj.trembleScore_csplus.Dur = arrayfun( @(i) (1/ (obj.ingress_index_csplus(i) - 2000 - obj.ingressBuffer) ).*sum(sum( any(im2bw( obj.cwt_csplus{i}( obj.cwt_window, [2000: obj.ingress_index_csplus(i) - obj.ingressBuffer ]), obj.threshold_min )), 1 )), [1:numel(obj.ingress_index_csplus)], 'ErrorHandler', @(x,y) NaN );
            obj.trembleScore_csminus.Dur = arrayfun( @(i) (1/ (obj.ingress_index_csminus(i) - 2000 - obj.ingressBuffer) ).*sum(sum( any(im2bw( obj.cwt_csminus{i}( obj.cwt_window, [2000: obj.ingress_index_csminus(i) - obj.ingressBuffer ]), obj.threshold_min )), 1 )), [1:numel(obj.ingress_index_csminus)], 'ErrorHandler', @(x,y) NaN );
            
            toc;
            
        end
        
        function getIngressScore( obj )
            fprintf('Getting ingress score for %s\n',obj.mouse);
            max2 = @(x) max(max(x));
            
            % This is the artificial number of trials for this mouse but
            % not the number of trials in the experiment overall, so it
            % will need to be padded at the end of this function 
            Ntrials = size( obj.displacement_csplus, 1 );
            
            aucs_ = arrayfun( @(trial) trapz( obj.displacement_csplus( trial , obj.ingress_index_csplus(trial):end ) ), [1:Ntrials] );
            
            % Currently the maximum AUC possible is calibrated to the
            % maximum displacement achieved by a given mouse
            aucs_rect = arrayfun( @(trial) [numel(obj.displacement_csplus(trial,:)) - obj.ingress_index_csplus(trial)] * max2(obj.displacement_csplus), [1:Ntrials], 'ErrorHandler', @(x,y) NaN );
            aucs_(aucs_<0) = 0;
            ingressScore_csplus = aucs_./aucs_rect;
            
            aucs_ = arrayfun( @(trial) trapz( obj.displacement_csminus( trial , obj.ingress_index_csminus(trial):end ) ), [1:Ntrials] );
            aucs_rect = arrayfun( @(trial) [numel(obj.displacement_csminus(trial,:)) - obj.ingress_index_csminus(trial)] * max2(obj.displacement_csplus), [1:Ntrials] );
            aucs_(aucs_<0) = 0;
            ingressScore_csminus = aucs_./aucs_rect;
            
            if numel( ingressScore_csplus )<obj.trials; ingressScore_csplus = [ ingressScore_csplus, nan(1, obj.trials - numel(ingressScore_csplus) ) ]; end
            if numel( ingressScore_csminus )<obj.trials; ingressScore_csminus = [ ingressScore_csminus, nan(1, obj.trials - numel(ingressScore_csminus) ) ]; end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compare the fit line to the ingress-model %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            corrcoefs_csplus = arrayfun( @(i) corrcoef( obj.fitlines_CSplus(i,:), obj.displacement_csplus(i,:) ), [1:obj.trials], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN );
            corrcoefs_csminus = arrayfun( @(i) corrcoef( obj.fitlines_CSminus(i,:), obj.displacement_csminus(i,:) ), [1:obj.trials], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN );
            
            corrcoefs_csplus = cellfun( @(i) i(2), corrcoefs_csplus, 'ErrorHandler', @(x,y) NaN );
            corrcoefs_csminus = cellfun( @(i) i(2), corrcoefs_csminus, 'ErrorHandler', @(x,y) NaN );
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add corrcoefs to the ingress score stats  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.ingressScore_csplus = table( corrcoefs_csplus', ingressScore_csplus', 'VariableNames', {'AUC_Fit','AUC_Mag'} );
            obj.ingressScore_csminus = table( corrcoefs_csminus', ingressScore_csminus', 'VariableNames', {'AUC_Fit','AUC_Mag'} );
            
        end
        
        function fitLine_allTrials( obj, varargin ) % Runs fitline on every trial, CSplus and CSminus
            
            fprintf('Fitting for Mouse %s\n', obj.mouse );
            trials = obj.trials; % Get all the trials
            if strcmp( varargin{1}, 'coeffs' ) % This will model the ingress only
                obj.ingress_stats_CSplus = arrayfun( @(x) obj.fitLine(x,'CSplus','coeffs'), [1:trials], 'UniformOutput', false )';
                obj.ingress_stats_CSminus = arrayfun( @(x) obj.fitLine(x,'CSminus','coeffs'), [1:trials], 'UniformOutput', false )';
            end
            
            if strcmp( varargin{1}, 'fitlines' ) % This will fit the whole trial and produce corrected timeseries data
                obj.fitlines_CSplus = cell2mat(arrayfun( @(x) obj.fitLine(x,'CSplus','fitlines'), [1:trials], 'UniformOutput', false )');
                obj.fitlines_CSminus = cell2mat(arrayfun( @(x) obj.fitLine(x,'CSminus','fitlines'), [1:trials], 'UniformOutput', false )');
                obj.corrected_displacement_csplus = obj.displacement_csplus - obj.fitlines_CSplus([1:size(obj.displacement_csplus,1)],:)
                obj.corrected_displacement_csminus = obj.displacement_csminus - obj.fitlines_CSminus([1:size(obj.displacement_csminus,1)],:);
            end
            
        end
        
        
        % Overloaded functions %
        function save( obj, varargin )
           
           
            if any(strcmp(varargin,'h5'))
                % Check if the directory exists
                myh5dir = sprintf('%s\\h5',cd);
                myh5file = sprintf('%s\\%s.h5',myh5dir,obj.mouse);
                
                if isdir( myh5dir )
                else
                   mkdir( myh5dir )
                end
                % Start saving
                
                if exist( myh5file )
                else
                    h5create(myh5file,'/displacement_csplus',size(obj.displacement_csplus) );
                    h5create(myh5file,'/displacement_csminus',size(obj.displacement_csminus) );
                    h5create(myh5file,'/ingress_index_csplus',size(obj.ingress_index_csplus) );
                    h5create(myh5file,'/ingress_index_csminus',size(obj.ingress_index_csminus) );
                end
                
                obj.writeH5displacements(myh5dir)
            end
            
        end
        
        function plot( obj, varargin )
           
            % Internal parameters for CWT %
            No = 10;
            Nv = 48;
            cwt = 0;
            ingress = 0;
            
            % Checking for options %
            if strcmp(varargin{1},'cwt')
                cwt = 1;
                w = arrayfun( @(trial) cwt( obj.corrected_displacement_csplus( trial , :), obj.fs, 'NumOctaves',No, 'VoicesPerOctave',Nv), [1:10], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN );
            else
                ingress = 1;
            end
            
            % Axes and tabs laid out for plotting all trials %
            if ingress
                figure('color','w'); 
                tabs_ = arrayfun( @(x) uitab('Title',sprintf('Trial %i',x)), [1:10] );
                ax_ = arrayfun( @(tab) axes(tab,'NextPlot','add','Ylim',[-3,20]), tabs_ );

                % Plotting %
                arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.displacement_csplus(x,:), 'b' ), [1:10] );
                arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.displacement_csminus(x,:), 'r' ), [1:10] );
                
                if ~isempty( obj.fitlines_CSplus )

                    arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.fitlines_CSplus(x,:), 'k' ), [1:10], 'ErrorHandler', @(x,y) [], 'UniformOutput', false );
                    %arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.fitlines_CSminus(x,:), 'k' ), [1:10], 'ErrorHandler', @(x,y) [], 'UniformOutput', false );

                end
                
            end
            
            % Axes and tabs laid out for plotting all trials %
            if cwt
                figure('color','w'); 
                tabs_ = arrayfun( @(x) uitab('Title',sprintf('Tab %i',x)), [1:10] );
                lims = size(w{1});
                im_ax_ = arrayfun( @(tab) axes(tab,'NextPlot','add','Ylim',[0,lims(1)],'Xlim',[0,lims(2)]), tabs_ );
                % Plotting %
                arrayfun(@(x) imagesc( w{x}.*conj(w{x}), 'parent', im_ax_(x),[0,0.003] ), [1:10] ); 
            end
                    
        end
       
        function varargout = changepts(varargin)

            window = 100;
            
            minpeakheight = varargin{ 1+find( strcmp( varargin, 'minpeakheight' )==1 ) };
            
            idx = [1:numel(varargin{2})];
            data = varargin{2};

                localz = zscore(arrayfun( @(x) mean(data(x:x+window))-mean(data(x-window:x)), [window+1:numel(data)-window] ));
                
                [pks,locs] = findpeaks(localz,'minpeakheight',prctile(localz,95),'minpeakdistance',100);
                varargout{1} = locs;
                varargout{2} = pks;
                
        end

    end
    
    methods(Static)
        
            function output = max_of_grad_adj_img( input_image )
           
                img = real( input_image );
                [~,Gy]= imgradientxy( img );
                %min_Gy = graythresh(Gy(:));
                %x = arrayfun( @(x) (x<min_Gy), Gy );
                featureIm = im2bw(imcomplement(Gy),1).*img;
                output = max(max(featureIm,[],1));
                if isempty(output); output = nan; end

            end
        
    end
    
end

