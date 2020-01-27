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
        dsFactor = 10;
    end
    
    properties( Access = private )
        t0 = 60000;
        tend = 200000;
        
        fourier_t0 = 2000;
        fourier_tend = 3000;
    end
    
    methods
        
        function obj = displacedMouse(mouseid, matfile_location)
            
            obj.mouse = mouseid;
            obj.file = matfile_location;
            
        end
        
        function associate_cs_index( obj, varargin )
            
            obj.csplus_third_index = varargin{1};
            obj.csminus_third_index = varargin{2};
            
        end
        
        function loadDisplacements( obj )
           
            object = load( obj.file );
            obj.displacement = object.displacement(:,[obj.t0:obj.dsFactor:obj.tend],:);
            
            if ~isempty(obj.csplus_third_index)
                obj.displacement_csplus = squeeze(obj.displacement(:,:,obj.csplus_third_index));
            end
            
            if ~isempty(obj.csminus_third_index)
                obj.displacement_csminus = squeeze(obj.displacement(:,:,obj.csminus_third_index));
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
            
            L = obj.fourier_tend - obj.fourier_t0;
            NFFT= 2^nextpow2( L );
            T = [obj.fourier_t0:obj.fourier_tend];
            fs=1000;
            fVals=fs*(0:NFFT/2-1)/NFFT;

            csplus_fft = arrayfun(@(x) fft(obj.displacement_csplus(x,T),NFFT), [1:size(obj.displacement_csplus,1)]', 'UniformOutput', false )
            csminus_fft = arrayfun(@(x) fft(obj.displacement_csminus(x,T),NFFT), [1:size(obj.displacement_csplus,1)]', 'UniformOutput', false )
            
            obj.csplus_fft_real = cellfun( @(X) X.*conj(X)/(NFFT*L), csplus_fft, 'UniformOutput', false ); %Power of each freq components
            obj.csminus_fft_real = cellfun( @(X) X.*conj(X)/(NFFT*L), csminus_fft, 'UniformOutput', false ); %Power of each freq components
            
            obj.tremblepower_csplus = cellfun( @(X) sum(X([5:12])), obj.csplus_fft_real ); % For a 1000 timepoint vector this is 4 Hz to 10 Hz
            obj.tremblepower_csminus = cellfun( @(X) sum(X([5:12])), obj.csminus_fft_real ); % For a 1000 timepoint vector this is 4 Hz to 10 Hz
            
        end
        
    end
end

