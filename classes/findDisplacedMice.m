classdef findDisplacedMice < handle
    % A class that works on all mouseObjs below
    
    properties
        folder
        keyTable
        ingressTable
        trembleSummary
        mouseObjs = struct();
        keyTableLocation
        ingressTableLocation 
        subfoldersTable
        mouseList
       ingress_score_summary_csplus
       ingress_score_summary_csminus
       tremble_score_summary_csplus
       tremble_score_summary_csminus
    end
    
    properties (Access=private)
       Ntrials = 10; 
    end
    
    % BEGINNING OF ORDINARY METHODS %
    
    methods
        function obj = findDisplacedMice(folder)
            % Looks for files in 'folder' that match the query P[A-Z](number)
            obj.folder = folder;
            
            folders = dir(folder);
            folders = folders( arrayfun( @(x) and( x.isdir, gt(numel(x.name),2)), folders ) );
            
            files = arrayfun( @(x) sprintf('%s\\%s\\DisplacementInOdor.mat',x.folder,x.name), folders, 'UniformOutput', false );
            
            allfiles = {};
            allids = {};
            
            for thisfile = files'
                mouseid = cell2mat( regexp( thisfile{1}, '[P][A-Z]\d+', 'match' ) );
                obj.mouseObjs.(mouseid) = displacedMouse( mouseid, thisfile{1} );
                allfiles = [allfiles; thisfile{1}];
                allids = [allids; mouseid];
            end
            
            obj.subfoldersTable = table( allfiles, allids );
            
        end
        
        function associateKey( obj, keyTableLocation )
% Looks into the provided 'keyTableLocation'. Populates the property keyTable, which contains information from cs_third_index_key.csv about which third index corresponds to CS+ and CS- 
            obj.keyTableLocation = keyTableLocation;
            obj.keyTable = readtable(keyTableLocation);
            
            for thisMouseObj = fields( obj.mouseObjs )'
                obj.mouseObjs.(thisMouseObj{1}).csplus_third_index = obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSp;
                obj.mouseObjs.(thisMouseObj{1}).csminus_third_index = obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSm;
            end
            
        end
         
        function fitLines_allMice( obj )
            structfun( @(x) x.fitLine_allTrials('trial'), obj.mouseObjs );
        end
        
        function fitIngress_allMice( obj )
            structfun( @(x) x.fitLine_allTrials('ingress'), obj.mouseObjs );
        end
        
        function getIngressScore_allMice( obj )
            structfun( @(x) x.getIngressScore(), obj.mouseObjs );
        end
        
        function getTrembleScore_allMice( obj )
            structfun( @(x) x.getTrembleScore(), obj.mouseObjs );
        end
        
        function associateIngressTimes( obj, ingressTableLocation )
% Looks into the provided 'ingressTableLocation'. Populates the property ingressTable, which contains ingress csplus frames (first ten columns) and csminus frames (second ten columns)  
            obj.ingressTableLocation = ingressTableLocation;
            
            obj.ingressTable = readtable( ingressTableLocation );
            
            for thisMouseObj = fields( obj.mouseObjs )'
                obj.mouseObjs.(thisMouseObj{1}).associate_csplus_index( table2array(obj.ingressTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.ingressTable.matched_names),...
                    find(cellfun(@(x) numel(regexp( x, 'csplus')), obj.ingressTable.Properties.VariableNames )==1)) ) );
                obj.mouseObjs.(thisMouseObj{1}).associate_csminus_index( table2array(obj.ingressTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.ingressTable.matched_names),...
                    find(cellfun(@(x) numel(regexp( x, 'csminus')), obj.ingressTable.Properties.VariableNames )==1)) ) );
            end
        end
        
        function loadDisplacements( obj ) 
% For each of mouseObjs in the object, this method runs loadDisplacements to extract data from the MAT file for each mouse, a method of the displacedMouse class.          
            for thisMouseObj = fields( obj.mouseObjs )'
                obj.mouseObjs.(thisMouseObj{1}).associate_cs_index( obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSp, ...
                                                                    obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSm );
            end
            
            structfun(@(x) x.loadDisplacements(), obj.mouseObjs );
        end
        
        
        function computeFourier( obj )
% For each of mouseObjs in the object, if there is loaded mouseObjs displacement data, this method uses data to run getFourier, a method of the displacedMouse class.  
            for thisMouseObj = fields( obj.mouseObjs )'
               obj.mouseObjs.(thisMouseObj{1}).getFourier(); 
            end
        end
        
        function makeTrembleSummary( obj )
% For each of mouseObjs in the object, if displacement has been analyzed using computeFourier, this method uses data to run computeFourier, a method of the displacedMouse class.  
            tmp_csplus = structfun( @(x) x.tremblepower_csplus, obj.mouseObjs, 'UniformOutput', false );
            tmp_csminus = structfun( @(x) x.tremblepower_csplus, obj.mouseObjs, 'UniformOutput', false );
            obj.trembleSummary.csplus = array2table(struct2array( tmp_csplus )',...
                'VariableNames',arrayfun(@(x) sprintf('Trial_%i',x), [1:obj.Ntrials],'UniformOutput',false),...
                'RowNames',fields(obj.mouseObjs));
            obj.trembleSummary.csminus = array2table(struct2array( tmp_csminus )',....
                'VariableNames',arrayfun(@(x) sprintf('Trial_%i',x), [1:obj.Ntrials],'UniformOutput',false),...
                'RowNames',fields(obj.mouseObjs));
        end
        
        function makeList( obj )
            obj.mouseList = fields( obj.mouseObjs );
        end
        
        function saveAsH5( obj, fieldsAsCell )
            for myfield = fieldsAsCell
                arrayfun( @(x) x.writeH5data( 'C:\\teamtremble_classes\\h5', myfield ), obj.mouseObjs );
            end
        end
        
        function makeScoreSummary( obj, varargin )
            
            if nargin>1
                trials = varargin{1};
            else
                trials = [1:obj.Ntrials];
            end
            
            tremble_score = structfun( @(x) x.trembleScore_csplus, obj.mouseObjs, 'UniformOutput', false );
            t_csplus_as_array = struct2array( tremble_score );
            obj.tremble_score_summary_csplus = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';
            obj.tremble_score_summary_csplus = log( obj.tremble_score_summary_csplus(:,trials) );
            
            tremble_score = structfun( @(x) x.trembleScore_csminus, obj.mouseObjs, 'UniformOutput', false );
            t_csminus_as_array = struct2array( tremble_score );
            obj.tremble_score_summary_csminus = reshape( t_csminus_as_array, numel(t_csminus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';
            obj.tremble_score_summary_csminus = log( obj.tremble_score_summary_csminus(:,trials) );
            
            ingress_score = structfun( @(x) x.ingressScore_csplus, obj.mouseObjs, 'UniformOutput', false );
            t_csplus_as_array = struct2array( ingress_score );
            obj.ingress_score_summary_csplus = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';
            obj.ingress_score_summary_csplus = obj.ingress_score_summary_csplus(:,trials);
            
            ingress_score = structfun( @(x) x.ingressScore_csminus, obj.mouseObjs, 'UniformOutput', false );
            t_csminus_as_array = struct2array( ingress_score );
            obj.ingress_score_summary_csminus = reshape( t_csminus_as_array, numel(t_csminus_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';
            obj.ingress_score_summary_csminus = obj.ingress_score_summary_csminus(:,trials);
            
        end
        
        function printScoreSummary( obj, I , T )
            
            ingress_cutoff = I;
            tremble_cutoff = T;
            fprintf('\nProbability of I and T | CS+: %1.2f\n',sum( and(ingress_score_summary_csplus>ingress_cutoff,tremble_score_summary_csplus>tremble_cutoff) ) )
            fprintf('Probability of I and T | CS-: %1.2f\n',sum( and(ingress_score_summary_csminus>ingress_cutoff,tremble_score_summary_csminus>tremble_cutoff) ) )
            fprintf('Probability of I | CS+: %1.2f\n',sum(ingress_score_summary_csplus>ingress_cutoff) )
            fprintf('Probability of I | CS-: %1.2f\n',sum(ingress_score_summary_csminus>ingress_cutoff) )
            fprintf('Probability of T | CS+: %1.2f\n',sum(tremble_score_summary_csplus>tremble_cutoff) )
            fprintf('Probability of T | CS-: %1.2f\n',sum(tremble_score_summary_csminus>tremble_cutoff) )
        
        end

        
        function boxplot( obj, trialType, scoreType )
            
            
            Nmice = numel( fields(obj.mouseObjs) );
            Lbls = fields( obj.mouseObjs );
            Ntrials = obj.Ntrials;
            
            figure('color','w','position',[270,480,920,240]);
            
            if strcmp(scoreType,'tremble')
                scoreType = sprintf('trembleScore_%s',trialType);
                tremble_score = structfun( @(x) x.(scoreType), obj.mouseObjs, 'UniformOutput', false );
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
                trembleType = sprintf('trembleScore_csplus');
                tremble_score = structfun( @(x) x.(trembleType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( tremble_score );
                tremble_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';
                cnts = histc( log( tremble_score(:) ),[0:1:6] );
                plot( [0:1:6], cnts/sum(cnts) );
                
                trembleType = sprintf('trembleScore_csminus');
                tremble_score = structfun( @(x) x.(trembleType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( tremble_score );
                tremble_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';
                cnts = histc( log( tremble_score(:) ),[0:1:6] );
                plot( [0:1:6], cnts/sum(cnts) );
                
                set(gca,'Clim',[0,6]);
                title('Tremble scores','FontWeight','Normal','FontSize',16);
                legend({'CS+','CS-'});
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
            
            if nargin<3
                fprintf('Too few inputs. Usage:\nplot( obj, trialType, scoreType )\n\nOptions for trialType are csplus, csminus\nOptions for scoreType are tremble or ingress');
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
                tremble_score = structfun( @(x) x.(trembleType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( tremble_score );
                tremble_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )';
                % Must flip up-down because of the labeling?
                imagesc( (log(tremble_score)), 'parent', ax );
                set(gca,'Clim',[0,6]);
                title('Tremble scores','FontWeight','Normal','FontSize',16);
            end
            
            if strcmp( scoreType, 'ingress' )
                ingressType = sprintf('ingressScore_%s',trialType);
                ingress_score = structfun( @(x) x.(ingressType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( ingress_score );
                ingress_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';
                % Must flip up-down because of the labeling?
                imagesc( ingress_score, 'parent', ax );
                set(gca,'Clim',[0,1]);
                title('Ingress scores','FontWeight','Normal','FontSize',16);
            end
            
            if strcmp( scoreType, 'both' )
                ingressType = sprintf('ingressScore_%s',trialType);
                ingress_score = structfun( @(x) x.(ingressType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( ingress_score );
                ingress_score = reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(ingress_score),1), size(fields(ingress_score),1) )';
                
                trembleType = sprintf('trembleScore_%s',trialType);
                tremble_score = structfun( @(x) x.(trembleType), obj.mouseObjs, 'UniformOutput', false );
                t_csplus_as_array = struct2array( tremble_score );
                tremble_score = log( reshape( t_csplus_as_array, numel(t_csplus_as_array)/size(fields(tremble_score),1), size(fields(tremble_score),1) )' );
                
                ax = ax
                Nmice = size(ingress_score, 1 );
                mycolors = lines( Nmice );
                arrayfun( @(x) plot( ingress_score(x,:), tremble_score(x,:), 'o', 'markerfacecolor', mycolors(x,:), 'markersize', 8, 'markeredgecolor', 'k' ), [1:Nmice] )
                
            end
            
            datacursormode on
            dcm_obj = datacursormode(fig);
            % NOTE THIS QUIRK: Had to use a static method for the dcm
            % function as ordinary methods didn't work
            set(dcm_obj, 'UpdateFcn', {@findDisplacedMice.myfxn,obj,trialType})
            
            if ~strcmp( scoreType, 'both' ); colorbar(); end
        end
        
    end
    % END OF ORDINARY METHODS %
    
    % BEGINNING OF STATIC METHODS %
    methods(Static)

        function txt = myfxn(~,event_obj,obj,trialType)
            trial = event_obj.Position(1);
            mouseNum = event_obj.Position(2); 
            if isempty( obj.mouseList ); obj.makeList(); end;
            figure('color','w'); 
            
            mouse_data = obj.mouseObjs.(obj.mouseList{mouseNum});
            
            ax = arrayfun( @(x) subplot(3,1,x,'Xlim',[0,14000],'NextPlot','add','Ylim',[-4,20]), [1:3] );
            if strcmp( trialType,'csplus' )
                plot( ax(1), mouse_data.displacement_csplus(trial,:) )
                line( [mouse_data.ingress_index_csplus(trial),mouse_data.ingress_index_csplus(trial)], [-4,20], 'color', 'k', 'parent', ax(1) );
                title( sprintf('%s Trial: %i (score=%1.2f)', obj.mouseList{mouseNum}, trial, mouse_data.trembleScore_csplus(trial)), 'parent', ax(1) );
                plot( ax(2), mouse_data.corrected_displacement_csplus(trial,:) )
                imagesc( mouse_data.cwt_csplus{trial}, 'parent', ax(3) )
                set(gca,'YLim',[0,481]);
                txt = mouse_data.trembleScore_csplus(trial);
            end
            
            if strcmp( trialType,'csminus' )
                plot( ax(1), mouse_data.displacement_csminus(trial,:) )
                line( [mouse_data.ingress_index_csminus(trial),mouse_data.ingress_index_csminus(trial)], [-4,20], 'color', 'k', 'parent', ax(1) );
                title( sprintf('%s Trial: %i (score=%1.2f)', obj.mouseList{mouseNum}, trial, mouse_data.trembleScore_csminus(trial)), 'parent', ax(1) );
                plot( ax(2), mouse_data.corrected_displacement_csminus(trial,:) )
                imagesc( mouse_data.cwt_csminus{trial}, 'parent', ax(3) )
                set(gca,'YLim',[0,481]);
                txt = mouse_data.trembleScore_csminus(trial);
            end
            
        end

    end
    % END OF STATIC METHODS %
        
end

