function labelTMApatients(curDataParams)
    
    [isOctave,matType,matExt] = checkIsOctave();

    folNamL = { 'TMA_JG010_010', ...
                'TMA_JG010_011', ...
                'TMA_JG010_012', ...
                'TMA_JG010_013', ...
                'TMA_JG010_014', ...
                'TMA_JG010_015'         };

    tmaNo = find(strcmp(curDataParams.folNam,folNamL));
    if isempty(tmaNo)
        disp('TMA patient labelling not supported')
    else
        
        dataTypeNam = dataTypeNamSelect('Binned',curDataParams);
        load(['./miscMATfiles/testerTMAdata_all' dataTypeNam],'p_number','tma')

        mFileNam = matFileNamSelect('Raw',curDataParams);
        load([mFileNam matExt],'LXY','fExists','emptySpec')
        fExists_all = fExists + emptySpec > 0;
        spec_raw = fExists(fExists_all);
        if curDataParams.smoothParam > 0
            mFileNam = matFileNamSelect('Binned',curDataParams);
            load([mFileNam matExt],'fExists')
        end
        spec = fExists(fExists_all);
        
        p_number = p_number(tma == tmaNo);
        temp = zeros(sum(spec_raw),1);
        temp(spec(spec_raw)) = p_number(spec_raw(spec));
        p_number = temp;

        p_list = unique(p_number);
        for p_idx = 1:length(p_list)
            p = p_list(p_idx);
            text(median(LXY(p_number==p,1)),median(LXY(p_number==p,2)),num2str(p));
        end
    end
        
    
    
end