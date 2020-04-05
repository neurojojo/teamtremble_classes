classdef displacedMouse < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mouse
        file
        displacement
        displacement_csplus
        displacement_csminus
        csplus_third_index
        csminus_third_index
        ingress_index_csplus
        ingress_index_csminus
        csplus_fft_real
        csminus_fft_real
        tremblepower_csplus
        tremblepower_csminus
        displacement_csplus_filtered
        displacement_csminus_filtered
        fitlines_CSplus
        fitlines_CSminus
        corrected_displacement_csplus
        corrected_displacement_csminus
        h5_datasets
        trembleScore_csplus
        trembleScore_csminus
        cwt_csplus
        cwt_csminus
        ingress_stats_CSplus
        ingress_stats_CSminus
        aucs_fraction
        ingressScore_csplus
        ingressScore_csminus
    end
    
    properties( Access = private )
        
        fs = 1000; % Downsampled sampling frequency
        dsFactor = 10;
        trials = 10;
        time_vector = (0:14000)/1000;
        
        Nframes = [];
        
        t0 = 60000;
        tend = 200000;
        
        d = 400; % An internal parameter that is calibrated to the duration before and after an ingress
        % Use this in the fitLine subfunction
        
        fourier_t0 = 2000;
        fourier_tend = 3000;
        
        % Butterworth filter properties
        butter_n = 4;
        butter_fpass = [(1/10000) 30];
        butter_a;
        butter_b;
        
    end
    
    methods
        
        function obj = displacedMouse(mouseid, matfile_location)
            
            obj.mouse = mouseid;
            obj.file = matfile_location;
            [obj.butter_b,obj.butter_a] = butter(obj.butter_n,2*obj.butter_fpass(1,:)*2/obj.fs);
            obj.Nframes = numel( obj.time_vector );
            
        end
        
        function associate_cs_index( obj, varargin )
            
            obj.csplus_third_index = varargin{1};
            obj.csminus_third_index = varargin{2};
            
        end
        
        function loadDisplacements( obj )
            fprintf('Loading displacement for %s\n', obj.mouse );
            
            object = load( obj.file );
            obj.displacement = object.displacement(:,[obj.t0:obj.dsFactor:obj.tend],:);
            
            if ~isempty(obj.csplus_third_index)
                obj.displacement_csplus = squeeze(obj.displacement(:,:,obj.csplus_third_index));
                fprintf('Filtering\n');
                obj.displacement_csplus_filtered = cell2mat( ...
                    arrayfun(@(x) filtfilt( obj.butter_b, obj.butter_a, obj.displacement_csplus(x,:) ), [1:obj.trials],...
                    'UniformOutput', false, 'ErrorHandler', @(x,y) zeros(1,size(obj.displacement_csplus,2)) )' );
            end
            
            if ~isempty(obj.csminus_third_index)
                obj.displacement_csminus = squeeze(obj.displacement(:,:,obj.csminus_third_index));
                obj.displacement_csminus_filtered = cell2mat( ...
                    arrayfun(@(x) filtfilt( obj.butter_b, obj.butter_a, obj.displacement_csminus(x,:) ), [1:obj.trials],...
                    'UniformOutput', false, 'ErrorHandler', @(x,y) zeros(1,size(obj.displacement_csminus,2)) )' );
            end
            
            if and(~isempty(obj.displacement_csminus),~isempty(obj.displacement_csplus))
                obj.displacement = [];
            end
            
        end
        
        function identifyIngress(obj)

            for i=1:size(for_figures,1)

                plot( for_figures(i,:).Displacement );
                lines_ = for_figures(i,:).Locs{1};

                rowfun( @(x,y) text(y,25,sprintf('%i',x) ), table( [1:numel(lines_)]', lines_ ) );
                arrayfun( @(x) line([x,x],[-20,30],'color','k'), lines_ );
                mychangept(i) = inputdlg('Test');

            end
            
        end
        
        function associate_csplus_index(obj, array_ingress_times)
            obj.ingress_index_csplus = array_ingress_times;
        end
        
        function associate_csminus_index(obj, array_ingress_times)
            obj.ingress_index_csminus = array_ingress_times;
        end
       
        function getFourier(obj)
% Computes csplus_fft_real, csminus_fft_real, tremblepower_csplus, tremblepower_csminus
            fprintf('Computing fourier for %s\n', obj.mouse );
            anyfilter = @(x) x;%-mean(x(:));
            
            L = obj.fourier_tend - obj.fourier_t0;
            NFFT= 2^nextpow2( L );
            T = [obj.fourier_t0:obj.fourier_tend];
            fs=obj.fs;
            fVals=fs*(0:NFFT/2-1)/NFFT;

            csplus_fft = arrayfun(@(x) fft( anyfilter(obj.displacement_csplus(x,T)),NFFT), [1:size(obj.displacement_csplus,1)]', 'UniformOutput', false );
            csminus_fft = arrayfun(@(x) fft( anyfilter(obj.displacement_csminus(x,T)),NFFT), [1:size(obj.displacement_csminus,1)]', 'UniformOutput', false );
            
            obj.csplus_fft_real = cellfun( @(X) X.*conj(X)/(NFFT*L), csplus_fft, 'UniformOutput', false ); %Power of each freq components
            obj.csminus_fft_real = cellfun( @(X) X.*conj(X)/(NFFT*L), csminus_fft, 'UniformOutput', false ); %Power of each freq components
            
            obj.tremblepower_csplus = cellfun( @(X) sum(X([5:12])), obj.csplus_fft_real ); % For a 1000 timepoint vector this is 4 Hz to 10 Hz
            obj.tremblepower_csminus = cellfun( @(X) sum(X([5:12])), obj.csminus_fft_real ); % For a 1000 timepoint vector this is 4 Hz to 10 Hz
            
        end

        % Line-fitting %
        function myfitline = fitLine( obj, trial, trialType, varargin ) % Fits a piecewise linear function
            d = 1000;
            Nframes = obj.Nframes;
            
            fprintf('Computing for trial %i\n', trial );
            
            % Check for whether the trialType is CS+ or CS- %
            if strcmp(trialType,'CS+')
                Y = obj.displacement_csplus(trial,:);
                ing = obj.ingress_index_csplus(trial);
                if eq(ing,0); myfitline = zeros(1,Nframes); return; end; % Create zero-vector for displacements lacking ingress
            else if strcmp(trialType,'CS-')
                Y = obj.displacement_csminus(trial,:);
                ing = obj.ingress_index_csminus(trial);
                if eq(ing,0); myfitline = zeros(1,Nframes); return; end; % Create zero-vector for displacements lacking ingress
                else
                    fprintf("Specify either 'CS+' or 'CS-' as a third argument\n");
                    myfitline = [];
                    return
                end
            end
               
            if strcmp( varargin{1}, 'ingress' )
                
                Y_ = Y;
                Y_([1:ing-d]) = 0;
                Y_([ing-d+1:min( numel(Y), ing+d )]) = Y([ing-d+1:min( numel(Y), ing+d )]);
                Y_(ing+d+1:end) = mean(Y(ing+d+1:end));

                [pks_min,locs_min] = min( Y_ );
                [pks_max,locs_max] = max( Y_ );

                % Parabolic fit
                coeffs = polyfit( [0:locs_max-locs_min], Y_(locs_min:locs_max), 2 );

                Fit_Table = struct();
                Fit_Table.X_pre = [0:locs_min-1];
                Fit_Table.Y_pre = zeros( numel(1, Fit_Table.X_pre ), 1 );
                Fit_Table.X_ing = [locs_min:locs_max-1];
                Fit_Table.Y_ing = polyval( coeffs, [locs_min:locs_max-1]-locs_min )';
                Fit_Table.X_post = [locs_max+1:numel(Y)];
                Fit_Table.Y_post = mean(Y([locs_max+1:numel(Y)])) * ones( numel(Fit_Table.X_post), 1 );

                y = [ Fit_Table.Y_pre; Fit_Table.Y_ing; Fit_Table.Y_post ];
                myfitline = struct();
                if locs_min > locs_max
                    fprintf('Error for %s on trial %i. Check ingress time!\n', obj.mouse, trial)
                    myfitline = nan;
                else
                    cc = corrcoef( y, Y ); cc = cc(2);
                    myfitline = cc;
                end
                
            end
            
            if strcmp( varargin{1}, 'trial' )
                myfitline = [ smooth( Y(1:ing), d ); Y(ing+1:end)' ]';
            end
            
        end
        
        function getTrembleScore( obj )
           
            fprintf('Computing tremble score for %s\n', obj.mouse );
            No = 10;
            Nv = 48;
            d = 800;
            
            % CS-PLUS
            if ~isempty(obj.cwt_csplus)
                w_conj = obj.cwt_csplus;
            else
                w = arrayfun( @(trial) cwt( obj.corrected_displacement_csplus( trial , :), obj.fs, 'NumOctaves',No, 'VoicesPerOctave',Nv), [1:10], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN );
                close all;
                w_conj = arrayfun( @(trial) w{trial}.*conj( w{trial} ), [1:10], 'UniformOutput', false );
            end
            
            pre_mean = arrayfun( @(trial) mean2(w_conj{trial}([200:350],[1000:2000])), [1:10], 'UniformOutput', true, 'ErrorHandler', @(x,y) NaN   );
            post_mean = arrayfun( @(trial) mean2(w_conj{trial}([200:350],[2000:obj.ingress_index_csplus(trial)-d])), [1:10], 'UniformOutput', true, 'ErrorHandler', @(x,y) NaN  );
            
            obj.cwt_csplus = w_conj;
            obj.trembleScore_csplus = post_mean./pre_mean;
            
            % CS-MINUS
            if ~isempty(obj.cwt_csminus)
                w_conj = obj.cwt_csminus;
            else
                w = arrayfun( @(trial) cwt( obj.corrected_displacement_csminus( trial , :), obj.fs, 'NumOctaves',No, 'VoicesPerOctave',Nv), [1:10], 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN );
                close all;
                w_conj = arrayfun( @(trial) w{trial}.*conj( w{trial} ), [1:10], 'UniformOutput', false );
            end
            
            pre_mean = arrayfun( @(trial) mean2(w_conj{trial}([200:350],[1000:2000])), [1:10], 'UniformOutput', true, 'ErrorHandler', @(x,y) NaN   );
            post_mean = arrayfun( @(trial) mean2(w_conj{trial}([200:350],[2000:obj.ingress_index_csplus(trial)-d])), [1:10], 'UniformOutput', true, 'ErrorHandler', @(x,y) NaN  );
            
            obj.cwt_csminus = w_conj;
            obj.trembleScore_csminus = post_mean./pre_mean;
            
        end
        
        function myinterpolatedline = interpLine( obj, trial, trialType ) % Fits a piecewise linear function
           
            if strcmp(trialType,'CS+')
                Y = obj.corrected_displacement_csplus(trial,:);
                ing = obj.ingress_index_csplus(trial);
            else if strcmp( trialType,'CS-');
            	Y = obj.corrected_displacement_csminus(trial,:);
                ing = obj.ingress_index_csplus(trial);
                end
            end
            
            D = abs( zscore( diff( Y ) ) );
            [~,sorted_D] = sort(D);
            interpolate_idx = [ sorted_D(end-1), sorted_D(end) ];
            interp_pts = numel([interpolate_idx(1):interpolate_idx(2)]);
            
            Y([interpolate_idx(1):interpolate_idx(2)]) = linspace( Y(interpolate_idx(1)), Y(interpolate_idx(2)), numel( Y([interpolate_idx(1):interpolate_idx(2)]) ) );
            %myfit = polyfit( [interpolate_idx(1):interpolate_idx(2)]-interpolate_idx(1), Y( [ interpolate_idx(1):interpolate_idx(2) ] ), 2 );
            
            %Y([interpolate_idx(1):interpolate_idx(2)]) = polyval( myfit, [interpolate_idx(1):interpolate_idx(2)]-interpolate_idx(1) );
            
            %interp1( Y([interpolate_idx(1)-100:interpolate_idx(1)-10, interpolate_idx(2)+10:interpolate_idx(2)+100]), ...
            %    [interpolate_idx(1)-100:interpolate_idx(1)-10, interpolate_idx(2)+10:interpolate_idx(2)+100], ...
            %   [interpolate_idx(1):interpolate_idx(2)] )
            
        end
        
        function getIngressScore( obj )
           
            max2 = @(x) max(max(x));
            
            aucs_ = arrayfun( @(trial) trapz( obj.displacement_csplus( trial , obj.ingress_index_csplus(trial):end ) ), [1:10] );
            %aucs_rect = arrayfun( @(trial) [numel(obj.displacement_csplus(trial,:)) - obj.ingress_index_csplus(trial)] * max(obj.displacement_csplus( trial , obj.ingress_index_csplus(trial):end )), [1:10] );
            aucs_rect = arrayfun( @(trial) [numel(obj.displacement_csplus(trial,:)) - obj.ingress_index_csplus(trial)] * max2(obj.displacement_csplus( : , obj.ingress_index_csplus(trial):end )), [1:10] );
            aucs_(aucs_<0) = 0;
            obj.ingressScore_csplus = aucs_./aucs_rect;
            
            aucs_ = arrayfun( @(trial) trapz( obj.displacement_csminus( trial , obj.ingress_index_csminus(trial):end ) ), [1:10] );
            %aucs_rect = arrayfun( @(trial) [numel(obj.displacement_csminus(trial,:)) - obj.ingress_index_csminus(trial)] * max(obj.displacement_csminus( trial , obj.ingress_index_csminus(trial):end )), [1:10] );
            aucs_rect = arrayfun( @(trial) [numel(obj.displacement_csplus(trial,:)) - obj.ingress_index_csplus(trial)] * max2(obj.displacement_csplus( : , obj.ingress_index_csplus(trial):end )), [1:10] );
            aucs_(aucs_<0) = 0;
            obj.ingressScore_csminus = aucs_./aucs_rect;
            
        end
        
        function fitLine_allTrials( obj, varargin ) % Runs fitline on every trial, CS+ and CS-
            
            fprintf('Fitting for Mouse %s\n', obj.mouse );
            trials = obj.trials; % Get all the trials
            if strcmp( varargin{1}, 'ingress' ) % This will model the ingress only
                obj.ingress_stats_CSplus = arrayfun( @(x) obj.fitLine(x,'CS+','ingress'), [1:trials], 'UniformOutput', true )';
                obj.ingress_stats_CSminus = arrayfun( @(x) obj.fitLine(x,'CS-','ingress'), [1:trials], 'UniformOutput', true )';
            end
            
            if strcmp( varargin{1}, 'trial' ) % This will fit the whole trial and produce corrected timeseries data
                obj.fitlines_CSplus = cell2mat(arrayfun( @(x) obj.fitLine(x,'CS+','trial'), [1:trials], 'UniformOutput', false )');
                obj.fitlines_CSminus = cell2mat(arrayfun( @(x) obj.fitLine(x,'CS-','trial'), [1:trials], 'UniformOutput', false )');
                obj.corrected_displacement_csplus = obj.displacement_csplus - obj.fitlines_CSplus;
                obj.corrected_displacement_csminus = obj.displacement_csminus - obj.fitlines_CSminus;
            end
            
        end
        
        function writeH5data( obj, myh5dir, field )
            
            myh5file = sprintf('%s\\%s.h5',myh5dir,obj.mouse);
            
            if exist( myh5file ) % h5 File already present
                myh5info = h5info( myh5file );
                datasetnames = arrayfun(@(x) x.Name, myh5info.Datasets , 'UniformOutput', false);
                if any(strcmp(datasetnames,field)) % Field already exists in the database
                    fprintf('The field %s already exists and is being overwritten\n',field);
                else
                    h5create( myh5file, sprintf('/%s',field), size( obj.(field) ) );
                end
            else % If file does not exist %
                if size( obj.(field) , 1 ) > 10
                    fprintf('Writing compressed h5\n');
                    h5create( myh5file, sprintf('/%s',field), size( obj.(field) ) );
                else
                    h5create( myh5file, sprintf('/%s',field), size( obj.(field) ) );
                end
            end
            
            h5write( myh5file, sprintf('/%s',field), obj.(field) );
            % Add the h5 info to the log structure
            obj.addH5toLog( myh5file, sprintf('/%s',field), datestr(datetime) );
            
        end
        
        function addH5toLog(obj, myh5file, myh5field, myh5)
            
            if ~istable( obj.h5_datasets ); obj.h5_datasets = table('VariableTypes',{'string','string','string'}, 'Size',[0,3],'VariableNames',{'H5_file','Field','Time_created'}); end;
            obj.h5_datasets(end+1,:) = [ { myh5file }, { myh5field }, { myh5 } ];
            
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
            
            % Checking for options %
            if strcmp(varargin{1},'cwt')
                    w = arrayfun( @(trial) cwt( obj.corrected_displacement_csplus( trial , :), obj.fs, 'NumOctaves',No, 'VoicesPerOctave',Nv), [1:10], 'UniformOutput', false );
            end
            close all;
            
            % Axes and tabs laid out for plotting all trials %
            figure; 
            tabs_ = arrayfun( @(x) uitab('Title',sprintf('Tab %i',x)), [1:10] );
            ax_ = arrayfun( @(tab) axes(tab,'NextPlot','add','Ylim',[-3,20]), tabs_ );
            
            
            % Plotting %
            arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.corrected_displacement_csplus(x,:) ), [1:10] );
            arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.corrected_displacement_csminus(x,:) ), [1:10] );
            
            if ~isempty( obj.fitlines_CSplus )
                
                arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.fitlines_CSplus(x,:), 'k' ), [1:10], 'ErrorHandler', @(x,y) [], 'UniformOutput', false );
                arrayfun(@(x) plot( ax_(x), obj.time_vector(1:end), obj.fitlines_CSminus(x,:), 'k' ), [1:10], 'ErrorHandler', @(x,y) [], 'UniformOutput', false );
                
            end
            
            
            % Axes and tabs laid out for plotting all trials %
            figure; 
            tabs_ = arrayfun( @(x) uitab('Title',sprintf('Tab %i',x)), [1:10] );
            lims = size(w{1});
            im_ax_ = arrayfun( @(tab) axes(tab,'NextPlot','add','Ylim',[0,lims(1)],'Xlim',[0,lims(2)]), tabs_ );
            % Plotting %
            arrayfun(@(x) imagesc( w{x}.*conj(w{x}), 'parent', im_ax_(x),[0,0.003] ), [1:10] );
                    
        end
        
    end
end

