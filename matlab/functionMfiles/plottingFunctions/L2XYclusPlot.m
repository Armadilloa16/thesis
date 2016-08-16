function map = L2XYclusPlot(y,fExists,swapXY,highlight)

    if length(y) ~= sum(sum(fExists))
        error('ERROR: !!! incorrect number of spectra input into L2XY')
    elseif sum(y==0)+sum(y==1)+sum(y==2)+sum(y==3)+sum(y==4)+sum(y==5)+sum(y==6) ~= sum(sum(fExists))
        error('ERROR: !!! incorrect format for 4clus')
    end
    
    colChoice = 1;
    
    switch colChoice
        case 1
            % 1 - Pink/Purple -- off tissue
            % 2 - Cyan -- adipose
            % 3 - Orange -- cancer
            % 4 - Grey/Purple.
            % 5 - Blue
            % 6 - Green -- stroma
            nCol = max(y);
            C = [100,80,80,10,80,80];
            L = [80,85,85,85,80,85];
            H = [300,180,60,270,240,120];
            C = C(1:nCol);
            L = L(1:nCol);
            H = H(1:nCol);
            if nargin > 3 
                if sum(y == highlight) > 0
                    L = L + 10;
                    C = C - 50;
                    L(highlight) = L(highlight) - 10;
                    C(highlight) = C(highlight) + 50; 
                else
                    error('no no no no no')
                end
            end
            RGB = hcl2rgb([H',C',L']);
    
        case 2
            % 1 - Pink/Purple -- off tissue
            % 2 - Cyan -- adipose
            % 3 - Orange -- cancer
            % 4 - Green -- stroma
            % 5 - Blue
            % 6 - Grey/Purple.

            RGB = [ 0.7 0.7 0.7;
                    1.0 0.0 0.0;
                    0.0 0.0 1.0;
                    0.3 0.3 1.0;
                    0.0 0.0 0.0;
                    1.0 0.0 1.0];
            nCol = max(y);
            RGB = RGB(1:nCol,:);    
            
        case 3
            
            RGB = [ 0.8 0.8 0.8;
                    1.0 0.0 0.0;
                    0.0 0.0 1.0;
                    0.3 1.0 0.3;
                    0.5 0.5 0.5;
                    1.0 1.0 1.0];
            nCol = max(y);
            RGB = RGB(1:nCol,:);    
            
        case 4
            
            nCol = max(y);
            
            RGB = repmat((0:((1-(1/nCol))/(nCol-1)):(1-(1/nCol)))',1,3);    
%             RGB = RGB([1 3 4 2],:);
            
    end
            
            
            
    yCols = ones(length(y),3);
    for idx = 1:nCol
        yCols(y==idx,:) = repmat(RGB(idx,:),sum(y==idx),1); 
    end
    
    map = ones(size(fExists,1),size(fExists,2),3);
    mat = ones(size(fExists,1),size(fExists,2));
    mat(fExists) = yCols(:,1);
    map(:,:,1) = mat;
    mat(fExists) = yCols(:,2);
    map(:,:,2) = mat;
    mat(fExists) = yCols(:,3);
    map(:,:,3) = mat;
  
    temp = map;
    if ~swapXY
        map = ones(size(fExists,2),size(fExists,1),3);
        for k = 1:3
            map(:,:,k) = temp(:,:,k)';
        end
    else
        map = temp;
    end
    
end