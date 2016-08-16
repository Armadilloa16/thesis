function [] = savePlot(plotParam,dataParam,formats_ind)
    
    outputFolnam = './Output/';
    formats = {'fig' 'pdf' 'eps' 'bmp' 'jpg' 'png'};
    genNoLabels = true;
    
    dataTypeNam = dataTypeNamSelect('Binned',dataParam);
    fnam = '';
                        
    switch plotParam.plotType
        case 'Clus'
            fnam = ['_' plotParam.distance num2str(plotParam.K)];
            if plotParam.highlight
                fnam = [fnam 'highlight' num2str(plotParam.k)];
            end
            
        case 'PCA'
            fnam = ['_PC' num2str(plotParam.K)];
            
        case 'Heatmap'
            fnam = ['_' plotParam.distance num2str(plotParam.K)];
            fnam = [fnam 'clus' num2str(plotParam.k)];
            if isfield(plotParam,'isLegend')
                if plotParam.isLegend
                    fnam = [fnam 'legend'];
                    genNoLabels = false;
                end
            end
                    
    end
    
    fnam = [plotParam.plotType '_' dataParam.folNam fnam dataTypeNam];
    if isfield(plotParam,'fixed')
        if plotParam.fixed
            fnam = [fnam '_fixed'];
        end
    end
    
    
    
    
    
    
    if nargin == 2
        formats_ind = [1 1 1 0 0 0] ~= 0;
    elseif nargin == 3
        if length(size(formats_ind)) == length(size(formats))
            if size(formats_ind,1) == size(formats,1) && ...
                    size(formats_ind,2) == size(formats,2)
                
                if iscell(formats_ind)
                    formats_ind = strcmp(formats_ind,formats);
                elseif isnumeric(formats_ind)
                    formats_ind = formats_ind ~= 0;
                elseif ~islogical(formats_ind)
                    error('format_ind is the wrong format')
                end
                
            else
                error(['formats_idx must be either a ' ...
                    num2str(size(formats,1)) ' x ' ...
                    num2str(size(formats,2)) ' cell or numeric ' ...
                    'indicator array'])
            end
        else
            error(['formats_idx must be either a ' ...
                    num2str(size(formats,1)) ' x ' ...
                    num2str(size(formats,2)) ' cell or numeric ' ...
                    'indicator array'])
        end
    else
        error('WAKAKAKAAKWHAAAAA?')
    end
        
    for format = formats(formats_ind)
        saveas(gcf,[outputFolnam fnam],format{1})
    end
    if genNoLabels
        title('')
        xlabel('')
        ylabel('')
        for format = formats(formats_ind)
            saveas(gcf,['./Output/' fnam '_noLabels'],format{1})
        end
    end
    
     
end