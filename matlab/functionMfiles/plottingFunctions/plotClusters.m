function plotClusters(clus,curDataParams,tnam,altCols,highlight)
    
    if nargin > 3
        if length(altCols) ~= max(clus)
            error('bad!')
        end
        temp = clus;
        for i = 1:length(altCols)
            clus(temp == i) = altCols(i);
        end
        clear temp
    end
       
    [~,~,matExt] = checkIsOctave();

    mFileNam = matFileNamSelect('Supplementary Data');
    load([mFileNam matExt],'metaData')
    invX = metaData.(curDataParams.folNam).invX;
    invY = metaData.(curDataParams.folNam).invY;
    swapXY = metaData.(curDataParams.folNam).swapXY;
    clear metaData
    
    mFileNam = matFileNamSelect('Raw',curDataParams);
    load([mFileNam matExt],'fExists','emptySpec','X','Y')
    if isfield(curDataParams,'smoothParam')
        if curDataParams.smoothParam > 0
            mFileNam = matFileNamSelect('Binned',curDataParams);
            load([mFileNam matExt],'fExists','emptySpec')
        end
    end
    clear emptySpec
        
    if nargin > 4
        clusMAP = L2XYclusPlot(clus,fExists,swapXY,highlight);
    else
        clusMAP = L2XYclusPlot(clus,fExists,swapXY);
    end
    
    if swapXY
        temp = invX;
        invX = invY;
        invY = temp;
        clear temp  
    end
    if ~swapXY
        image(X,Y,clusMAP)
        squarePlot(X,Y)
        xlabel('X')
        ylabel('Y')
    else
        image(Y,X,clusMAP)
        squarePlot(Y,X)
        xlabel('Y')
        ylabel('X')
    end
    
    if nargin > 3
        title(tnam);
    end
    if invX
        set(gca,'Xdir','reverse')
    end
    if invY
        set(gca,'Ydir','normal')
    end

end