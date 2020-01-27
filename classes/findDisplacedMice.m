classdef findDisplacedMice < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        folder
        keyTable
        ingressTable
        trembleSummary
        mouseObjs = struct();
        keyTableLocation
        ingressTableLocation 
    end
    
    properties (Access=private)
       Ntrials = 10; 
    end
    
    methods
        function obj = findDisplacedMice(folder)
            
            obj.folder = folder;
            
            folders = dir(folder);
            folders = folders( arrayfun( @(x) and( x.isdir, gt(numel(x.name),2)), folders ) );
            
            files = arrayfun( @(x) sprintf('%s\\%s\\DisplacementInOdor.mat',x.folder,x.name), folders, 'UniformOutput', false );
            
            for thisfile = files'
                mouseid = cell2mat( regexp( thisfile{1}, '[P][A-Z]\d+', 'match' ) );
                obj.mouseObjs.(mouseid) = displacedMouse( mouseid, thisfile{1} );
            end
            
        end
        
        function associateKey( obj, keyTableLocation ) % A three column CSV
            obj.keyTableLocation = keyTableLocation;
            obj.keyTable = readtable(keyTableLocation);
            
            for thisMouseObj = fields( obj.mouseObjs )'
                obj.mouseObjs.(thisMouseObj{1}).csplus_third_index = obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSp;
                obj.mouseObjs.(thisMouseObj{1}).csminus_third_index = obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSm;
            end
        end
         
        function associateIngressTimes( obj, ingressTableLocation )
            obj.ingressTableLocation = ingressTableLocation;
            
            obj.ingressTable = readtable( ingressTableLocation );
            
            for thisMouseObj = fields( obj.mouseObjs )'
                obj.mouseObjs.(thisMouseObj{1}).associate_csplus_index( table2array(obj.ingressTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.ingressTable.matched_names),...
                    find(cellfun(@(x) numel(regexp( x, 'csplus')), obj.ingressTable.Properties.VariableNames )==1)) ) );
                obj.mouseObjs.(thisMouseObj{1}).associate_csminus_index( table2array(obj.ingressTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.ingressTable.matched_names),...
                    find(cellfun(@(x) numel(regexp( x, 'csminus')), obj.ingressTable.Properties.VariableNames )==1)) ) );
            end
        end
        
        function loadDisplacements( obj ) % A three column CSV
            
            for thisMouseObj = fields( obj.mouseObjs )'
                obj.mouseObjs.(thisMouseObj{1}).associate_cs_index( obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSp, ...
                                                                    obj.keyTable( cellfun(@(x) strcmp(x,thisMouseObj), obj.keyTable.Animal), :).CSm );
            end
            
            structfun(@(x) x.loadDisplacements(), obj.mouseObjs );
        end
        
        
        function computeFourier( obj )
            for thisMouseObj = fields( obj.mouseObjs )'
               obj.mouseObjs.(thisMouseObj{1}).getFourier(); 
            end
        end
        
        function makeTrembleSummary( obj )
            tmp_csplus = structfun( @(x) x.tremblepower_csplus, obj.mouseObjs, 'UniformOutput', false );
            tmp_csminus = structfun( @(x) x.tremblepower_csplus, obj.mouseObjs, 'UniformOutput', false );
            obj.trembleSummary.csplus = array2table(struct2array( tmp_csplus )',...
                'VariableNames',arrayfun(@(x) sprintf('Trial_%i',x), [1:obj.Ntrials],'UniformOutput',false),...
                'RowNames',fields(obj.mouseObjs));
            obj.trembleSummary.csminus = array2table(struct2array( tmp_csminus )',....
                'VariableNames',arrayfun(@(x) sprintf('Trial_%i',x), [1:obj.Ntrials],'UniformOutput',false),...
                'RowNames',fields(obj.mouseObjs));
        end
        
    end
end

