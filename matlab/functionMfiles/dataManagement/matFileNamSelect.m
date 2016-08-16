function [mFileNam,regexpVars] = matFileNamSelect(dataType,paramStruct)
    
%     EXAMPLE USE:
% 
%     [mFileNam,regexpVars] = matFileNamSelect('Clus',curDataParams);
%     save(matType,[mFileNam matExt],'-regexp',regexpVars)
% 

    matFileFolder = './dataMATfiles/';

    switch dataType
        case 'Dataset Index'
            mFileNam = [matFileFolder 'datasetNames'];
            regexpVars = 'dataNams|dataTypeNams|dataExists';
            
        case 'Supplementary Data'
            mFileNam = [matFileFolder 'supplementaryData'];
            regexpVars = 'useAltSeed4clus|metaData';
            
        case 'Etma Clinical Data'
            mFileNam = [matFileFolder 'Etma_clinical_data'];
            regexpVars = 'p_list|p_lnm|p_suit|p_lvsi|p_grade|p_size|p_FIGO_old|p_FIGO_old_lvls|p_FIGO_new|p_FIGO_new_lvls|p_T_old|p_T_old_lvls|p_T_new|p_T_new_lvls|p_myo_inv|p_myo_thi|p_ser_dis';
            
        case 'Vtma Clinical Data'
            mFileNam = [matFileFolder 'Vtma_clinical_data'];
            regexpVars = 'p_list|p_lnm|p_suit|p_grade';

        otherwise
            dataTypeNam = dataTypeNamSelect(dataType,paramStruct);
            mFileNam = [matFileFolder paramStruct.folNam dataTypeNam];
            
            switch dataType
                case 'Raw'
                    regexpVars = 'L|X|Y|R|emptySpec|fExists|LXY|XYL|Vars';

                case 'Annotation'
                    regexpVars = 'p_number|annot_cancer';

                case 'Binned'
                    if ~isfield(paramStruct,'smoothParam')
                        if ~isfield(paramStruct,'dataType')
                            regexpVars = 'mcdata|vbincentrs|mdataL';
                        elseif strcmp(paramStruct.dataType,'Binary')
                            regexpVars = 'mcdata|vbincentrs|mdataL';
                        else
                            regexpVars = 'mdata';
                        end
                    elseif paramStruct.smoothParam == 0
                        if ~isfield(paramStruct,'dataType')
                            regexpVars = 'mcdata|vbincentrs|mdataL';
                        elseif strcmp(paramStruct.dataType,'Binary')
                            regexpVars = 'mcdata|vbincentrs|mdataL';
                        else
                            regexpVars = 'mdata';
                        end
                    else
                        regexpVars = 'mbdata|vbincentrs|fExists|emptySpec|nIter|converged';                
                    end
                    
                case 'Normalisation'
                    regexpVars = 'Shat|nCal|meanCalIntensity|CV';

                case {'PCA','pPCA'}
                    regexpVars = 'veigval|meigvec|vmean|mpc';
                    
                case 'Clus'
                    if ~isfield(paramStruct,'dataType')
                        regexpVars = 'clus|centroicell|DIPScell|cutOFFcell|plotCutOFFcell|plotDistcell|Kmax';
                    elseif strcmp(paramStruct.dataType,'Binary')
                        regexpVars = 'clus|centroicell|DIPScell|cutOFFcell|plotCutOFFcell|plotDistcell|Kmax';
                    else
                        regexpVars = 'clus|centroicell|Kmax';
                    end
                    
                case 'rSummary'
                    if ~isfield(paramStruct,'dataType')
                        regexpVars = 'rdata|nSpec|r_list|r_list_count';
                    elseif strcmp(paramStruct.dataType,'Binary')
                        regexpVars = 'rdata|nSpec|r_list|r_list_count';
                    else
                        regexpVars = 'rdata';
                    end
                        
                case 'pSummary'
                    if ~isfield(paramStruct,'dataType')
                        regexpVars = 'pdata|p_list|p_lnm|p_suit|nSpec|vbincentrs';
                    elseif strcmp(paramStruct.dataType,'Binary')
                        regexpVars = 'pdata|p_list|p_lnm|p_suit|nSpec|vbincentrs';
                    else
                        regexpVars = 'pdata';
                    end
                    
                case 'pClassification'
                    regexpVars = 'class_out|dvec|b|class_out_LOO|dvec_LOO|b_LOO';
                    
            end
           
    end 
    
end