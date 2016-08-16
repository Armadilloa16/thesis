function dataTypeNam = dataTypeNamSelect(dataType,paramStruct)

    switch dataType
        case 'Raw'
            dataTypeNam = '_rawData';
                                    
        case 'Annotation'
            dataTypeNam = '_annotation';

        case 'Normalisation'
            if paramStruct.ppmNormTol
                dataTypeNam = ['_tol' num2str(paramStruct.normTol) 'ppm'];
            else
                dataTypeNam = ['_tol' num2str(1000*paramStruct.normTol) 'mDa'];
            end    
            switch paramStruct.dataType
                case {'Intensity','LogIntensity'}
                    dataTypeNam = [dataTypeNam '_intensity'];
                case 'Area'
                    dataTypeNam = [dataTypeNam '_area'];
                case 'SN'
                    dataTypeNam = [dataTypeNam '_sn'];
            end
            dataTypeNam = [dataTypeNam '_normalisation'];
            
        case {'Binned','PCA','pPCA','Clus','Intensity','Area','SN','rSummary','pSummary','pClassification'}
            dataTypeNam = ['_Bin' num2str(100*paramStruct.binSize)];
            if isfield(paramStruct,'wiggle')
                if paramStruct.wiggle
                    if isfield(paramStruct,'wiggle_den')
                        if paramStruct.wiggle_den == 2
                            dataTypeNam = [dataTypeNam '_wiggle'];                            
                        elseif isfield(paramStruct,'wiggle_num')
                            if paramStruct.wiggle_num ~= 0
                                dataTypeNam = [dataTypeNam '_' num2str(paramStruct.wiggle_num) 'wiggle' num2str(paramStruct.wiggle_den)];
                            end
                        end
                    else
                        if ~isfield(paramStruct,'wiggle_num')
                            dataTypeNam = [dataTypeNam '_wiggle'];
                        else
                            error('wiggle_den exists but wiggle_num does not!')
                        end
                    end
                end
            end
            if isfield(paramStruct,'smoothParam')
                if paramStruct.smoothParam ~= 0
                    dataTypeNam = [dataTypeNam '_' num2str(100*paramStruct.smoothParam) 'smooth'];
                end
            end
            if strcmp(dataType,'rSummary')
                if strcmp(paramStruct.dataType,'LogIntensity')
                    dataTypeNam = [dataTypeNam '_log'];
                end
            end
            if strcmp(dataType,'pSummary')
                if strcmp(paramStruct.dataType,'LogIntensity')
                    dataTypeNam = [dataTypeNam '_log'];
                end
            end
            if strcmp(dataType,'pPCA')
                if strcmp(paramStruct.dataType,'LogIntensity')
                    dataTypeNam = [dataTypeNam '_log'];
                end
            end
            if strcmp(dataType,'pClassification')
                if strcmp(paramStruct.dataType,'LogIntensity')
                    dataTypeNam = [dataTypeNam '_log'];
                end
            end
            if isfield(paramStruct,'dataType')
                switch paramStruct.dataType
                    case {'LogIntensity' 'Intensity'}
                        dataTypeNam = [dataTypeNam '_intensity'];
                    case 'Area'
                        dataTypeNam = [dataTypeNam '_area'];
                    case 'SN'
                        dataTypeNam = [dataTypeNam '_sn'];
                end
            end         
            switch dataType                    
                case 'PCA'
                    dataTypeNam = [dataTypeNam '_pca'];
                case 'Clus'
                    dataTypeNam = [dataTypeNam '_' paramStruct.clusType 'Clus'];
                case {'rSummary','pSummary','pClassification','pPCA'}
                    if paramStruct.normalisation
                        dataTypeNam = [dataTypeNam '_normalised'];
                    end
                    dataTypeNam = [dataTypeNam '_minNcal' num2str(paramStruct.nCal_tol)];
                    if paramStruct.use_cancer_annot
                        dataTypeNam = [dataTypeNam '_cancerAnnotatedSpectraOnly'];
                    end
                    switch dataType
                        case 'pSummary'
                            dataTypeNam = [dataTypeNam '_minNspecPerRegion' num2str(paramStruct.nSpec_r_tol)];
                            dataTypeNam = [dataTypeNam '_pSummary'];
                        case 'rSummary'
                            dataTypeNam = [dataTypeNam '_rSummary'];
                        case 'pClassification'
                            dataTypeNam = [dataTypeNam '_minNspecPerRegion' num2str(paramStruct.nSpec_r_tol)];
                            dataTypeNam = [dataTypeNam '_minNspecPerBin' num2str(paramStruct.nSpec_bin_tol)];
                            if ~paramStruct.includeEmptyVals
                                dataTypeNam = [dataTypeNam '_noEmptyVals'];
                            end                            
                            dataTypeNam = [dataTypeNam '_' paramStruct.classificationMethod];
                            switch paramStruct.classificationMethod
                                case {'pcaNB' 'pcaLDA' 'pcaDWD'}
                                    dataTypeNam = [dataTypeNam '_' num2str(paramStruct.nComponents) 'pcs'];
                                case {'cca1NB' 'cca1LDA' 'cca1DWD' 'cca2NB' 'cca2LDA' 'cca2DWD'}
                                    dataTypeNam = [dataTypeNam '_' num2str(paramStruct.nComponents) 'var'];
                            end
                            dataTypeNam = [dataTypeNam '_pClassification'];
                        case 'pPCA'
                            dataTypeNam = [dataTypeNam '_minNspecPerRegion' num2str(paramStruct.nSpec_r_tol)];
                            dataTypeNam = [dataTypeNam '_minNspecPerBin' num2str(paramStruct.nSpec_bin_tol)];
                            if ~paramStruct.includeEmptyVals
                                dataTypeNam = [dataTypeNam '_noEmptyVals'];
                            end
                            if paramStruct.restrict_p_suit
                                dataTypeNam = [dataTypeNam '_onlySuitP'];
                            end
                            dataTypeNam = [dataTypeNam '_pPCA'];
                    end
            end
            
        otherwise
            error('mFileNam Selection went awry!')
    
    end

end