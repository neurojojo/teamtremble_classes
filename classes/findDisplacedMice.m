classdef findDisplacedMice < handle
    % A class that works on all mouseObjs below
    
    properties
        folder
        keyTableLocation
        keyTable
        ingressTable
        trembleSummary
        mouseObjs = struct()
        ingressTableLocation 
        subfoldersTable
        mouseList
        ingress_score_summary_csplus
        ingress_score_summary_csminus
        tremble_score_summary_csplus
        tremble_score_summary_csminus
        stimorder_struct
        
        all_stim_orders
        all_stim_orders_str
        
    end
    
    properties (Access=private)
       Ntrials = 10; % Number of trials per cs type 
       CStypes = 2;
       Nmice;
    end
    
    % BEGINNING OF ORDINARY METHODS %
    
    methods
        function obj = findDisplacedMice(folder)
            % Looks for files in 'folder' that match the query P[A-Z](number)
            obj.folder = folder;
            
            folders = dir(folder);
            folders = folders( arrayfun( @(x) and( x.isdir, gt(numel(x.name),2)), folders ) );
            
            folders = arrayfun( @(x) sprintf('%s\\%s\\',x.folder,x.name), folders, 'UniformOutput', false );
            
            tableVars = {'mouseId','displacement_file','displacement_exist','stimulus_file','stimulus_exist'};
            
            obj.subfoldersTable = table( [],[],[],[],[], 'VariableNames',tableVars );
            
            for thisfolder = folders'
                mouseid = cell2mat( regexp( thisfolder{1}, '[P][A-Z]\d+', 'match' ) );
                if ~isempty(mouseid)
                    obj.mouseObjs.(mouseid) = displacedMouse( mouseid, thisfolder{1} );
                    stim_file = fullfile( thisfolder{1}, 'stimuli.mat' );
                    displacement_file = fullfile( thisfolder{1}, 'DisplacementInOdor.mat' );
                    obj.subfoldersTable = [ obj.subfoldersTable; ...
                        cell2table( {mouseid, displacement_file, exist(displacement_file), stim_file, exist(stim_file) },'VariableNames',tableVars) ];
                end
            end
            
            
            
        end
        
        function selectIngress_allMice( obj )
             structfun( @(x) x.selectIngress(), obj.mouseObjs );
        end
        
        function associateKey( obj, keyTableLocation )
            
            %keyTableLocation = fullfile(obj.datadir,'stimuli.mat');
            if eq(0,exist(keyTableLocation)); 
                fprintf('The stimulus key file is not located where you specified'); return; 
            end
            
            obj.keyTable = readtable(keyTableLocation);
            obj.keyTable.Animal = categorical( obj.keyTable.Animal );
            
            for thisMouseObj = fields( obj.mouseObjs )'
                
                % Stimuli.mat contains a array of integers corresponding to the identity of the stimuli
                stimulus_order_file = fullfile( obj.mouseObjs.(thisMouseObj{1}).datadir, 'stimuli.mat' );
                stimuli = load( stimulus_order_file );
                stimuli = stimuli.stimuli;
                %if ~strcmp(class( stimuli ),'single'); fprintf('Stimuli must be an array of numbers!'); return; end;
                
                obj.mouseObjs.(thisMouseObj{1}).csplus_index = [obj.keyTable( obj.keyTable.Animal == thisMouseObj, : ).CSp; ...
                                                                obj.keyTable( obj.keyTable.Animal == thisMouseObj, : ).CSp2];
                                                            
                % Last column of keyTable is a string in MATLAB 2018 so there needs to be an alternative version or catch the error
                obj.mouseObjs.(thisMouseObj{1}).csminus_index = [obj.keyTable( obj.keyTable.Animal == thisMouseObj, : ).CSm; ...
                                                                obj.keyTable( obj.keyTable.Animal == thisMouseObj, : ).CSm2];
                                                            
                
                obj.mouseObjs.(thisMouseObj{1}).keyTable = obj.keyTable( obj.keyTable.Animal == thisMouseObj, : );
                
                % All this does is use keyTable to decode the 1,2,3,4 from the stimulus file, where these numbers correspond to CSp,CSp2,CSm,CSm2 or O3 
                csp1 = table2array(obj.mouseObjs.(thisMouseObj{1}).keyTable(1,'CSp'));
                csp2 = table2array(obj.mouseObjs.(thisMouseObj{1}).keyTable(1,'CSp2'));
                csm1 = table2array(obj.mouseObjs.(thisMouseObj{1}).keyTable(1,'CSm'));
                csm2 = table2array(obj.mouseObjs.(thisMouseObj{1}).keyTable(1,'CSm2'));
                
                Nstimuli = max( stimuli );
                Ntrials = size( stimuli, 1 );
                
                switch Nstimuli
                    case 4
                        translation_table = table( [csp1;csp2;csm1;csm2],[{'p1'};{'p2'};{'m1'};{'m2'}], 'VariableNames', {'keyCode','stimType'} );
                    case 3
                        translation_table = table( [csp1;csm1;3],[{'p1'};{'m1'};{'o3'}], 'VariableNames', {'keyCode','stimType'} );
                end
                
                translation_table = sortrows( translation_table, 'keyCode' );
                
                empty_cell = cell( Ntrials , 1 );
                
                for keys_ = [1:4]
                    [~,locs] = find( stimuli==keys_ );
                    for i = find(stimuli==keys_)'
                        empty_cell{i} = translation_table(keys_,:).stimType{1};
                    end
                end
                
                obj.mouseObjs.(thisMouseObj{1}).keyString = sprintf('%s',categorical(empty_cell));
                
                for types_ = translation_table.stimType'
                    fprintf('Analyzing %s\n',thisMouseObj{1});
                    obj.mouseObjs.(thisMouseObj{1}).trialStruct.(types_{1}) = zeros(numel(stimuli),1);
                    thisType = 1+ ( regexp( obj.mouseObjs.(thisMouseObj{1}).keyString, types_{1} )-1 )/2;
                    obj.mouseObjs.(thisMouseObj{1}).trialStruct.(types_{1})( thisType ) = [ 1 : obj.Ntrials]';
                end
                rownames = arrayfun( @(x) sprintf('Trial_%i',x), [1:Ntrials], 'UniformOutput', false )
                obj.mouseObjs.(thisMouseObj{1}).trialTable = struct2table( obj.mouseObjs.(thisMouseObj{1}).trialStruct, 'RowNames', rownames );
                
                % Two experimental set-ups (1-3)
                if numel(translation_table.keyCode)==3
                   tmp1 = obj.mouseObjs.(thisMouseObj{1}).trialTable.p1;
                   tmp2 = obj.mouseObjs.(thisMouseObj{1}).trialTable.m1;
                   obj.mouseObjs.(thisMouseObj{1}).trialTable.p(tmp1~=0) = [ 1:obj.Ntrials ];
                   obj.mouseObjs.(thisMouseObj{1}).trialTable.m(tmp2~=0) = [ 1:obj.Ntrials ];
                end
                
                % Alternative experimental set-up (1-4)
                if numel(translation_table.keyCode)==4
                   tmp1 = obj.mouseObjs.(thisMouseObj{1}).trialTable.p1+obj.mouseObjs.(thisMouseObj{1}).trialTable.p2;
                   tmp2 = obj.mouseObjs.(thisMouseObj{1}).trialTable.m1+obj.mouseObjs.(thisMouseObj{1}).trialTable.m2;
                   obj.mouseObjs.(thisMouseObj{1}).trialTable.p(tmp1~=0) = [1:20];
                   obj.mouseObjs.(thisMouseObj{1}).trialTable.m(tmp2~=0) = [1:20];
                end
                
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Functions that are invoked only from a demo %
        % script or command line                      %
        %                                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fitLines_allMice( obj, type ) % Type can be 'coeffs' or 'fitlines'
            structfun( @(x) x.fitLine_allTrials(type), obj.mouseObjs );
        end
        
        function getIngressScore_allMice( obj )
            structfun( @(x) x.getIngressScore(), obj.mouseObjs );
        end
        
        function getTrembleScore_allMice( obj )
            structfun( @(x) x.getTrembleScore(), obj.mouseObjs );
        end
        
        function associateIngressTimes( obj, ingressTableLocation )
% Looks into the provided 'ingressTableLocation'. Populates the property ingressTable, which contains ingress csplus frames (first ten columns) and csminus frames (second ten columns)  
            if exist('ingressTableLocation')
                obj.ingressTableLocation = ingressTableLocation;
                obj.ingressTable = readtable( ingressTableLocation );

                for thisMouseObj = fields( obj.mouseObjs )'
                    mouserow = find( cellfun(@(x) strcmp(x,thisMouseObj), obj.ingressTable.mouse) == 1 );
                    tmp1 = table2array(obj.ingressTable( mouserow,...
                        find(cellfun(@(x) numel(regexp( x, 'csplus')), obj.ingressTable.Properties.VariableNames )==1)) ); tmp1=tmp1(isnan(tmp1)==0);
                    tmp2 = table2array(obj.ingressTable( mouserow,...
                        find(cellfun(@(x) numel(regexp( x, 'csminus')), obj.ingressTable.Properties.VariableNames )==1)) ); tmp2=tmp2(isnan(tmp2)==0);
                    obj.mouseObjs.(thisMouseObj{1}).associate_csplus_index( tmp1 );
                    obj.mouseObjs.(thisMouseObj{1}).associate_csminus_index(tmp2 );

                    obj.mouseObjs.(thisMouseObj{1}).ingress_index_csplus( isnan(obj.mouseObjs.(thisMouseObj{1}).ingress_index_csplus)==1 ) = size(obj.mouseObjs.(thisMouseObj{1}).displacement_csplus,2);
                    obj.mouseObjs.(thisMouseObj{1}).ingress_index_csminus( isnan(obj.mouseObjs.(thisMouseObj{1}).ingress_index_csminus)==1 ) = size(obj.mouseObjs.(thisMouseObj{1}).displacement_csminus,2);
                end
            else
                for thisMouseObj = fields( obj.mouseObjs )'
                    obj.mouseObjs.(thisMouseObj{1}).ingressTable.Type = categorical(obj.mouseObjs.(thisMouseObj{1}).ingressTable.Type);
                    obj.mouseObjs.(thisMouseObj{1}).ingress_index_csplus = obj.mouseObjs.(thisMouseObj{1}).ingressTable.Changepoint( obj.mouseObjs.(thisMouseObj{1}).ingressTable.Type=='csplus' );
                    obj.mouseObjs.(thisMouseObj{1}).ingress_index_csminus = obj.mouseObjs.(thisMouseObj{1}).ingressTable.Changepoint( obj.mouseObjs.(thisMouseObj{1}).ingressTable.Type=='csminus' );
                end
            end
            
        end
        
        function loadDisplacements_allMice( obj ) 
% For each of mouseObjs in the object, this method runs loadDisplacements to extract data from the MAT file for each mouse, a method of the displacedMouse class.                     
            structfun(@(x) x.loadDisplacements(), obj.mouseObjs );
            obj.mouseList = fields( obj.mouseObjs );
        end
        
        function result_table = search( obj, query )
            
            % Strip numbers %
            plus_or_minus = query( regexp(query,'[pm]') );
            len_sequence = numel( plus_or_minus ); % Tells you how many items are in a matching pattern %
            
            result = structfun( @(x) searchMouse(x,query), obj.mouseObjs , 'UniformOutput', false);
            fields_ = fields(result);
            empty_ = find( structfun( @(x) isempty(x), result )==1 );
            for field_ = fields_(empty_)
                result = rmfield(result,field_);
            end
            
            % Rearrange result for each mouse
            global_table = table();
            
            for field_ = fields(result)'
                mymat = table2array( result.(field_{1}) )';
                [numinseq,occurance] = find(mymat~=0);
                trialNum = mymat(:);
                mouseId = repmat( field_, numel(trialNum), 1 );
                trialType = repmat( plus_or_minus', max(occurance), 1 );
                tmp_table = table( mouseId,trialNum,numinseq,occurance,trialType );
                global_table = [global_table; tmp_table];
            end
            
            % After mice with no matching queries are discarded %
            tmpstruct = struct();
            result_table = table();
            for i = 1:size( global_table, 1 ); % Iterate over rows of tmp_table
                
                trial = global_table(i,:).trialNum;
                mouseId = global_table(i,:).mouseId{1};
                
                switch global_table(i,:).trialType
                    case 'p'
                        tmpstruct.AUC_Fit = obj.mouseObjs.(mouseId).ingressScore_csplus(trial,:).AUC_Fit;
                        tmpstruct.AUC_Mag = obj.mouseObjs.(mouseId).ingressScore_csplus(trial,:).AUC_Mag;
                        tmpstruct.Tremble_Mag = obj.mouseObjs.(mouseId).trembleScore_csplus.Mag(trial);
                        tmpstruct.Tremble_Dur = obj.mouseObjs.(mouseId).trembleScore_csplus.Dur(trial);

                    case 'm'
                        tmpstruct.AUC_Fit = obj.mouseObjs.(mouseId).ingressScore_csminus(trial,:).AUC_Fit;
                        tmpstruct.AUC_Mag = obj.mouseObjs.(mouseId).ingressScore_csminus(trial,:).AUC_Mag;
                        tmpstruct.Tremble_Mag = obj.mouseObjs.(mouseId).trembleScore_csminus.Mag(trial);
                        tmpstruct.Tremble_Dur = obj.mouseObjs.(mouseId).trembleScore_csminus.Dur(trial);
                end
                                
                        result_table = [result_table; struct2table( tmpstruct ) ];
                        
            end
            
            result_table = [global_table,result_table];
            
        end
        
        
        function mytable = summarize( obj )
            
            %results.ingress_mean = [ structfun( @(x) mean(x.ingressScore_csplus), obj.mouseObjs );...
            %    structfun( @(x) mean(x.ingressScore_csminus), obj.mouseObjs ) ];
            
            %results.tremble_mean = [ structfun( @(x) mean(x.trembleScore_csplus), obj.mouseObjs );...
            %    structfun( @(x) mean(x.trembleScore_csminus), obj.mouseObjs ) ];
            
            % Num of trials for each mouse
            N_csplus_trials_per_mouse = structfun( @(x) size(x.displacement_csplus,1), obj.mouseObjs );
            N_csminus_trials_per_mouse = structfun( @(x) size(x.displacement_csminus,1), obj.mouseObjs );
            Mice = fields(obj.mouseObjs);
            
            % Creating an array for the mouse names
            mouseIds = []; for i = 1:numel(Mice); mouseIds = [ mouseIds; repmat( {Mice{i}}, N_csplus_trials_per_mouse(i), 1 ) ]; end 
            for i = 1:numel(Mice); mouseIds = [ mouseIds; repmat( {Mice{i}}, N_csminus_trials_per_mouse(i), 1 ) ]; end 
            
            % Creating an array for the trial types
            Ntrials_Per_Type = 10;
            TrialType = []; for i = 1:numel(Mice); TrialType = [ TrialType; repmat( {'CS+1'}, Ntrials_Per_Type, 1 ); repmat( {'CS+2'}, N_csplus_trials_per_mouse(i)-Ntrials_Per_Type, 1 ) ]; end 
            for i = 1:numel(Mice); TrialType = [ TrialType; repmat( {'CS-1'}, Ntrials_Per_Type, 1 ); repmat( {'CS-2'}, N_csminus_trials_per_mouse(i)-Ntrials_Per_Type, 1 ) ]; end 
            
            TrialType_Generic = []; for i = 1:numel(Mice); TrialType_Generic = [ TrialType_Generic; repmat( {'CS+'}, N_csplus_trials_per_mouse(i), 1 ) ]; end 
            for i = 1:numel(Mice); TrialType_Generic = [ TrialType_Generic; repmat( {'CS-'}, N_csminus_trials_per_mouse(i), 1 ) ]; end 
            
            % Create an array for the trial number
            TrialNumber = [repmat( [1:Ntrials_Per_Type]', sum(N_csplus_trials_per_mouse)/Ntrials_Per_Type, 1 );...
                repmat( [1:Ntrials_Per_Type]', sum(N_csminus_trials_per_mouse)/Ntrials_Per_Type, 1 )];
            
            % Collect ingress magnitude %
            ingress = table( [flat(struct2array( structfun( @(x) x.ingressScore_csplus.AUC_Mag', obj.mouseObjs , 'UniformOutput', false)));...
                             flat(struct2array( structfun( @(x) x.ingressScore_csminus.AUC_Mag', obj.mouseObjs , 'UniformOutput', false)))] , 'VariableNames', {'Ingress_Mag'} );
            
            % Collect ingress onset %
            onset = table( [flat(struct2array(structfun( @(x) x.ingress_fraction_csplus, obj.mouseObjs , 'UniformOutput', false )));...
                flat(struct2array(structfun( @(x) x.ingress_fraction_csminus, obj.mouseObjs , 'UniformOutput', false )))] , 'VariableNames', {'Ingress_onset'} );
            
            % Collect tremble magnitude %
            tremble_mag = table( [flat(struct2array( structfun( @(x) x.trembleScore_csplus.Mag, obj.mouseObjs , 'UniformOutput', false)));...
                             flat(struct2array( structfun( @(x) x.trembleScore_csminus.Mag, obj.mouseObjs , 'UniformOutput', false)))] , 'VariableNames', {'Tremble_Mag'} );
            
            % Collect tremble duration %
            tremble_dur = table( [flat(struct2array( structfun( @(x) x.trembleScore_csplus.Dur, obj.mouseObjs , 'UniformOutput', false)));...
                             flat(struct2array( structfun( @(x) x.trembleScore_csminus.Dur, obj.mouseObjs , 'UniformOutput', false)))] , 'VariableNames', {'Tremble_Dur'} );
            
            % Collect tremble absolute duration %
            tremble_absdur = table( [flat(struct2array( structfun( @(x) x.trembleScore_csplus.AbsDur, obj.mouseObjs , 'UniformOutput', false)));...
                             flat(struct2array( structfun( @(x) x.trembleScore_csminus.AbsDur, obj.mouseObjs , 'UniformOutput', false)))] , 'VariableNames', {'Tremble_AbsDur'} );
            
            
            
            %fields_ = [cell2table( reshape( repmat( fields( obj.mouseObjs ), 1, obj.Ntrials.*2 )', 2*obj.Ntrials*obj.Nmice, 1 ), 'VariableNames', {'Mice'} );...
            %    cell2table( reshape( repmat( fields( obj.mouseObjs ), 1, obj.Ntrials.*2 )', 2*obj.Ntrials*obj.Nmice, 1 ), 'VariableNames', {'Mice'} )];
            
            %trials = array2table( repmat( [1:obj.Ntrials]', numel(fields_)/obj.Ntrials , 1 ), 'VariableNames', {'Trials'} );
            
            %trialtype = array2table( [repmat( [ repmat( {'CS+1'}, obj.Ntrials, 1 ); repmat( {'CS+2'}, obj.Ntrials, 1 ) ], Nmice, 1 );...
            %    repmat( [ repmat( {'CS-1'}, obj.Ntrials, 1 ); repmat( {'CS-2'}, obj.Ntrials, 1 ) ], Nmice, 1 )], 'VariableNames', {'TrialType'} );
            
            %trialtype_generic = categorical( [repmat( [ repmat( {'CS+'}, obj.Ntrials, 1 ); repmat( {'CS+'}, obj.Ntrials, 1 ) ], Nmice, 1 );...
            %    repmat( [ repmat( {'CS-'}, obj.Ntrials, 1 ); repmat( {'CS-'}, obj.Ntrials, 1 ) ], Nmice, 1 )] );
            
            %trialtype = array2table( repmat(  [repmat( {'CS+1'}, 10, 1 );...
            %                        repmat( {'CS+2'}, 10, 1 );...
            %                        repmat( {'CS-1'}, 10, 1 );...
            %                        repmat( {'CS-2'}, 10, 1 )] , 20, 1 ) , 'VariableNames', {'TrialType'} );
                                
            mytable = [ table(mouseIds), table(TrialType), table(TrialType_Generic), table(TrialNumber),...
                ingress,...
                tremble_mag,...
                tremble_dur,...
                tremble_absdur,...
                onset];
            mytable.TrialType = categorical( mytable.TrialType );
            mytable.TrialType_Generic = categorical( mytable.TrialType_Generic );
            mytable.mouseIds = categorical( mytable.mouseIds );
            
        end
        
        function boxplot( obj, trialType, scoreType )
            
            Nmice = numel( fields(obj.mouseObjs) );
            Lbls = fields( obj.mouseObjs );
            Ntrials = obj.Ntrials;
            
            figure('color','w','position',[270,480,920,240]);
            
            if strcmp(scoreType,'tremble')
                scoreType = sprintf('trembleScore_%s',trialType);
                tremble_score = structfun( @(x) x.(scoreType).Dur, obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( tremble_score );
                tremble_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';

                ax=subplot(1,5,1);
                boxplot( log(tremble_score) );
                set(gca,'XTick',[1:1:Ntrials],'XTickLabel',[0:1:Ntrials],'XLim',[0,Ntrials+1],...
                'NextPlot','add','box','off','YTick',[-6:1:6],'YTickLabel',[-6:1:6],'YAxisLocation','left',...
                'YTickLabelRotation',45,'YLim',[-6,6],'TickDir','out');
                xlabel('Trial number');
                box off;
                title( sprintf('Tremble (%s) by trial',trialType), 'FontWeight', 'normal'  );

                ax=subplot(1,5,[2:3]); 
                boxplot( log(tremble_score(:,[1:5]))');
                set(gca,'XTick',[1:1:Nmice],'XTickLabel',obj.mouseList,'XTickLabelRotation',60,...
                    'XLim',[0,Nmice+1],'NextPlot','add','box','off','YTick',[-6:1:6],'YTickLabel',[-6:1:6],'YAxisLocation','left',...
                'YTickLabelRotation',45,'YLim',[-6,6],'TickDir','out')
                box off;
                line([0,20],[1,1],'color','k')
                title( sprintf('Tremble (%s) by animal (first 5 trials)',trialType), 'FontWeight', 'normal' );
                
                ax=subplot(1,5,[4:5]); 
                boxplot( log(tremble_score(:,[6:10]))');
                set(gca,'XTick',[1:1:Nmice],'XTickLabel',obj.mouseList,'XTickLabelRotation',60,...
                    'XLim',[0,Nmice+1],'NextPlot','add','box','off','YTick',[-6:1:6],'YTickLabel',[-6:1:6],'YAxisLocation','left',...
                'YTickLabelRotation',45,'YLim',[-6,6],'TickDir','out')
                box off;
                line([0,20],[1,1],'color','k')
                title( sprintf('Tremble (%s) by animal (last 5 trials)',trialType), 'FontWeight', 'normal' );
            end
            
            if strcmp(scoreType,'ingress')
                
                scoreType = sprintf('ingressScore_%s',trialType);
                ingress_score = structfun( @(x) x.(scoreType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( ingress_score );
                ingress_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';

                ax=subplot(1,3,1);
                boxplot( ingress_score );
                set(gca,'XTick',[1:1:Ntrials],'XTickLabel',[1:1:Ntrials],'XLim',[0,Ntrials+1],...
                'NextPlot','add','box','off','YTick',[0:.1:1],'YTickLabel',[0:.1:1],'YAxisLocation','left',...
                'YTickLabelRotation',45,'YLim',[0,1],'TickDir','out');
                xlabel('Trial number');
                box off;
                title( sprintf('Ingress (%s) by trial',trialType), 'FontWeight', 'normal'  );

                ax=subplot(1,3,[2:3]); 
                boxplot( ingress_score' );
                set(gca,'XTick',[1:1:Nmice],'XTickLabel',obj.mouseList,'XTickLabelRotation',60,...
                    'XLim',[0,Nmice+1],'NextPlot','add','box','off','YTick',[0:.1:1],'YTickLabel',[0:.1:1],'YAxisLocation','left',...
                'YTickLabelRotation',45,'YLim',[0,1],'TickDir','out')
                box off;
                title( sprintf('Ingress (%s) by animal',trialType), 'FontWeight', 'normal' );
            end
            
        end
        
        function hist( obj, trialType, scoreType)
            
           if nargin<3
                fprintf('Too few inputs. Usage:\nhist( obj, trialType, scoreType )\n\nOptions for trialType are csplus, csminus\nOptions for scoreType are tremble or ingress\n');
                return
            end
            
            Nmice = numel( fields(obj.mouseObjs) );
            Lbls = fields( obj.mouseObjs );
            Ntrials = obj.Ntrials;
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Lay out the figures %
            %%%%%%%%%%%%%%%%%%%%%%%
            xlbl = arrayfun(@(x) sprintf('Trial %i',x), [1:Ntrials], 'UniformOutput', false );
            fig = figure('color','w'); 
            
            ax=axes('XTick',[0:1:6],'XTickLabel',[0:1:6],'XLim',[0,6],...
            'NextPlot','add','box','off','YTick',[0:.1:1],'YTickLabel',[0:.1:1],'YAxisLocation','left',...
            'YTickLabelRotation',45,'YLim',[0,1]);
            
            if strcmp( scoreType, 'tremble' )
                
                figure('color','w');
                set(gcf,'windowstate','maximized');
                drawnow();
                
                ax = arrayfun(@(x) subplot(2,4,x,'XColor','w','YColor','w','NextPlot','add','Units','pixels'),[1:8]);
                
                % CS-plus
                trembleType = sprintf('trembleScore_csplus');
                tremble_score_by_animal = structfun( @(x) x.(trembleType), obj.mouseObjs, 'UniformOutput', false );
                durations = struct2array( structfun(@(x) x.Dur, tremble_score_by_animal , 'UniformOutput', false) )';
                magnitudes = struct2array( structfun(@(x) x.Mag, tremble_score_by_animal , 'UniformOutput', false) )';
                
                % Subplot containing the scatter plot %
                scatter( durations, log(magnitudes),'parent',ax(5),'filled','markerfacealpha',0.7 );
                set( ax(5), 'XColor','k','YColor','k', 'YLim',[-10,5]);
                
                % Subplot to the right %
                tmp_ = histc( log(magnitudes), linspace(-10,5,20) ); tmp_ = tmp_./numel(magnitudes);
                bar( ax(6), linspace(-10,5,20), tmp_  ); 
                set( ax(6),'XLim',[-10,5],'YLim',[0,0.15],'XDir','reverse'); 
                camroll(ax(6),-90);
                set( ax(6), 'OuterPosition', get(ax(6),'OuterPosition')-[ 50 0 150 0] )
                height_ = get(ax(6),'OuterPosition') 
                height_ = height_(3);
                
                % Subplot to the top %
                tmp_ = histc( durations, linspace(0,1,20) ); tmp_ = tmp_./numel(durations);
                bar( ax(1), linspace(0,1,20), tmp_  ); 
                set( ax(1),'XLim',[-0.05,1],'YLim',[0,0.15] );
                pos_ = get(ax(1),'OuterPosition');
                %set( ax(1), 'OuterPosition', [pos_(1) pos_(2) pos_(3) height_] );
                set( ax(1),'OuterPosition', get(ax(1),'OuterPosition') - [0 90 0 0] )
                
                % CS-minus
                trembleType = sprintf('trembleScore_csminus');
                tremble_score_by_animal = structfun( @(x) x.(trembleType), obj.mouseObjs, 'UniformOutput', false );
                durations = struct2array( structfun(@(x) x.Dur, tremble_score_by_animal , 'UniformOutput', false) )';
                magnitudes = struct2array( structfun(@(x) x.Mag, tremble_score_by_animal , 'UniformOutput', false) )';
                
                % Subplot containing the scatter plot %
                scatter( durations, log(magnitudes),'parent',ax(7),'filled','markerfacealpha',0.7 );
                set( ax(7), 'XColor','k','YColor','k', 'YLim',[-10,5]);
                
                % Subplot to the right %
                tmp_ = histc( log(magnitudes), linspace(-10,5,20) ); tmp_ = tmp_./numel(magnitudes);
                bar( ax(8), linspace(-10,5,20), tmp_  ); 
                set( ax(8),'XLim',[-10,5],'YLim',[0,0.15],'XDir','reverse'); 
                camroll(ax(8),-90);
                set( ax(8), 'OuterPosition', get(ax(8),'OuterPosition')-[ 50 0 150 0] )
                height_ = get(ax(8),'OuterPosition')
                height_ = height_(3);
                
                % Subplot to the top %
                tmp_ = histc( durations, linspace(0,1,20) ); tmp_ = tmp_./numel(durations);
                bar( ax(3), linspace(0,1,20), tmp_  ); 
                set( ax(3),'XLim',[-0.05,1],'YLim',[0,0.15] );
                pos_ = get(ax(3),'OuterPosition');
                %set( ax(3), 'OuterPosition', [pos_(1) pos_(2) pos_(3) height_] );
                set( ax(3),'OuterPosition', get(ax(3),'OuterPosition') - [0 90 0 0] )
                
                ax(1).OuterPosition(4) = 250;
                ax(3).OuterPosition(4) = 250;
                
                set(ax(5),'ylabel',text(0,0,'Max tremble power'),'xlabel',text(0,0,'Fraction of time spent trembling') );
                set(ax(7),'ylabel',text(0,0,'Max tremble power'),'xlabel',text(0,0,'Fraction of time spent trembling') );
                
            end
            
            if strcmp( scoreType, 'ingress' )
                ingressType = sprintf('ingressScore_%s',trialType);
                ingress_score = structfun( @(x) x.(ingressType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( ingress_score );
                ingress_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';
                % Must flip up-down because of the labeling?
                hist( ingress_score, 'parent', ax );
                set(gca,'Clim',[0,1]);
                title('Ingress scores','FontWeight','Normal','FontSize',16);
            end
            
            
        end
        
        function plot( obj, trialType, scoreType )
            
            trialType = lower(trialType);
            
            if nargin<3
                fprintf('\nToo few inputs. Usage:\nplot( obj, trialType, scoreType )\n\nOptions for trialType are csplus, csminus\nOptions for scoreType are tremble or ingress');
                return
            end
            
            Nmice = numel( fields(obj.mouseObjs) );
            Lbls = fields( obj.mouseObjs );
            Ntrials = obj.Ntrials;
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Lay out the figures %
            %%%%%%%%%%%%%%%%%%%%%%%
            xlbl = arrayfun(@(x) sprintf('Trial %i',x), [1:Ntrials], 'UniformOutput', false );
            fig = figure('color','w'); 
            
            if ~strcmp( scoreType, 'both' )
                ax=axes('YTick',[1:1:Nmice],'YTickLabel',Lbls,'YLim',[.5,Nmice+.5],...
                'NextPlot','add','box','off','XTick',[1:1:Ntrials],'XTickLabel',xlbl,'XAxisLocation','top',...
                'XTickLabelRotation',45,'XLim',[.5,.5+Ntrials]); 
            else
                ax=axes('XTick',[0:0.1:1],'XTickLabel',[0:0.1:1],'XLim',[0,1],...
                'NextPlot','add','box','off','YTick',[-6:1:6],'YTickLabel',[-6:1:6],'YAxisLocation','left',...
                'YTickLabelRotation',45,'YLim',[-6,6]); 
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Decide what to plot (trembleScore or ingressScore) %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp( scoreType, 'tremble' )
                trembleType = sprintf('trembleScore_%s',trialType);
                tremble_score = structfun( @(x) x.(trembleType).Dur, obj.mouseObjs, 'UniformOutput', false );
                t_as_array = struct2array( tremble_score );
                tremble_score = reshape( t_as_array, numel(t_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';
                % Must flip up-down because of the labeling?
                imagesc( (tremble_score), 'parent', ax );
                set(gca,'Clim',[0,1]);
                title('Tremble durations','FontWeight','Normal','FontSize',16);
            end
            
            if strcmp( scoreType, 'ingress' )
                ingressType = sprintf('ingressScore_%s',trialType);
                ingress_score = structfun( @(x) x.(ingressType).AUC_Mag, obj.mouseObjs, 'UniformOutput', false );
                t_as_array = struct2array( ingress_score );
                ingress_score = reshape( t_as_array, numel(t_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';
                % Must flip up-down because of the labeling?
                imagesc( ingress_score, 'parent', ax );
                set(gca,'Clim',[0,1]);
                title('Ingress scores','FontWeight','Normal','FontSize',16);
            end
            
            datacursormode on
            dcm_obj = datacursormode(fig);
            % NOTE THIS QUIRK: Had to use a static method for the dcm
            % function as ordinary methods didn't work
            set(dcm_obj, 'UpdateFcn', {@findDisplacedMice.dynamic_plot_fxn,obj,trialType})
            
            if ~strcmp( scoreType, 'both' ); colorbar(); end
        end
        
    end
    % END OF ORDINARY METHODS %
    
    % BEGINNING OF STATIC METHODS %
    methods(Static)
        
        % Called by plot %
        function txt = dynamic_plot_fxn(~,event_obj,obj,trialType)
            trial = event_obj.Position(1);
            mouseNum = event_obj.Position(2); 
            figure('color','w'); 
            
            mouse_data = obj.mouseObjs.(obj.mouseList{mouseNum});
            
            ax = arrayfun( @(x) subplot(2,1,x,'Xlim',[0,14000],'NextPlot','add','Ylim',[-4,20],'TickDir','out'), [1:2] );
            if strcmp( trialType,'csplus' )
                plot( ax(1), mouse_data.displacement_csplus(trial,:) )
                line( [mouse_data.ingress_index_csplus(trial),mouse_data.ingress_index_csplus(trial)], [-4,20], 'color', 'k', 'parent', ax(1) );
                title( sprintf('%s Trial: %i (score=%1.2f)', obj.mouseList{mouseNum}, trial, mouse_data.trembleScore_csplus.Mag(trial)), 'parent', ax(1) );
                %plot( ax(2), mouse_data.corrected_displacement_csplus(trial,:) )
                imagesc( log(mouse_data.cwt_csplus{trial}), 'parent', ax(2), [-6,-1] )
                set(gca,'YLim',[0, size( mouse_data.cwt_csplus{trial},1 ) ]);
                txt = sprintf( 'Tremble magnitude: %1.2f',mouse_data.trembleScore_csplus.Mag(trial) );
            end
            
            if strcmp( trialType,'csminus' )
                plot( ax(1), mouse_data.displacement_csminus(trial,:) )
                line( [mouse_data.ingress_index_csminus(trial),mouse_data.ingress_index_csminus(trial)], [-4,20], 'color', 'k', 'parent', ax(1) );
                title( sprintf('%s Trial: %i (score=%1.2f)', obj.mouseList{mouseNum}, trial, mouse_data.trembleScore_csminus.Mag(trial)), 'parent', ax(1) );
                plot( ax(2), mouse_data.corrected_displacement_csminus(trial,:) )
                imagesc( mouse_data.cwt_csminus{trial}, 'parent', ax(3) )
                set(gca,'YLim',[0,481]);
                txt = sprintf( 'Tremble magnitude: %1.2f',mouse_data.trembleScore_csminus.Mag(trial) );
            end
            
            arrayfun(@(x) set(x,'XLabel',text(0,0,'Frame')), ax )
            set( ax(1), 'YLabel', text(0,0,'Displacement (mm)') );
            set( ax(2), 'YLabel', text(0,0,'Frequency'), 'YTick', [] );
            
        end
        
    end
    % END OF STATIC METHODS %
        
end

