function [class_out,d,b] = maldiDA(X,class,method,Xtest)
    
    % X         : a dxn data matrix
    % class     : a length n vector of class indices
    % method    : DA method to be used
    % Xtext     : dxm data matrix for testing (optional)
    %               -- if absent will default to Xtest = X
    
    
%     disp([num2str(sum(sum(X(:,class>0),2)==0)) ' empty variables removed.'])
%     disp([num2str(sum(sum(X(:,class>0),2)==sum(class>0))) ' full variables removed.'])
%     disp(' ')
    b = nan;
    valid_vars = var(X(:,class>0),0,2)>eps;
    if sum(valid_vars) == 0
        d = nan(size(valid_vars));
        if nargin == 4
            class_out = nan(1,size(Xtest,2));
            return
        end        
        class_out = nan(1,size(X,2));
        return
    end
    if nargin == 4
        Xtest = Xtest(valid_vars,:);
    end
    X = X(valid_vars,:);
    
    n_class = max(class);

    % Train
    switch method

        case 'LDA'
            Xbar = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                Xbar(:,class_idx) = mean(X(:,class==class_idx),2);
            end
            Xbarbar = mean(Xbar,2);
            Bi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                Bi(:,:,class_idx) = (Xbar(:,class_idx) - Xbarbar)*((Xbar(:,class_idx) - Xbarbar)');
            end
            B = sum(Bi,3);
            clear Bi
            Wi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                n = sum(class == class_idx);
                Wi(:,:,class_idx) = (1/(n-1))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))'));
            end
            W = sum(Wi,3);
            clear Wi
            if rcond(W) > eps
                [d,~] = eigs(W\B,1);
            else
                d = nan(size(W,1),1);
            end

        case 'NaiveBayes'            
            Xbar = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                Xbar(:,class_idx) = mean(X(:,class==class_idx),2);
            end
            Xbarbar = mean(Xbar,2);
            Bi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                Bi(:,:,class_idx) = (Xbar(:,class_idx) - Xbarbar)*((Xbar(:,class_idx) - Xbarbar)');
            end
            B = sum(Bi,3);
            clear Bi
            Wi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                n = sum(class == class_idx);
                Wi(:,:,class_idx) = diag((1/(n-1))*...
                    sum(((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n)).*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n)))),2));
            end
            W = sum(Wi,3);
            clear Wi
            [d,~] = eigs(W\B,1);

        case 'medoidLDA_Hamming'
            X = double(X);
            n = size(X,2);
            % Cosine Distance
    %         D = X'*X;
    %         D = D./repmat(sqrt(diag(D)),1,n);
    %         D = D./repmat(sqrt(diag(D))',n,1);
            % Euclidean Distance
    %         v = dot(X,X,1);
    %         D = bsxfun(@plus,v,v')-2*(X'*X);
            % Hamming Distance
            D = X'*X + (1-X)'*(1-X);
            D = (size(X,1) - D)./size(X,1);
            Xbar = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                [~,idx] = min(sum(D(:,class==class_idx),2));
                Xbar(:,class_idx) = double(X(:,idx));
            end
            Xbarbar = mean(Xbar,2);
            Bi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                Bi(:,:,class_idx) = (Xbar(:,class_idx) - Xbarbar)*((Xbar(:,class_idx) - Xbarbar)');
            end
            B = sum(Bi,3);
            clear Bi
            Wi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                n = sum(class == class_idx);
                Wi(:,:,class_idx) = (1/(n-1))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))'));
            end
            W = sum(Wi,3);
            clear Wi
            [V,D] = eig(W\B);
            [~,I] = max(abs(real(diag(D))));
            d = V(:,I);

        case 'medoidBayes_Hamming'
            X = double(X);
            n = size(X,2);
            % Cosine Distance
    %         D = X'*X;
    %         D = D./repmat(sqrt(diag(D)),1,n);
    %         D = D./repmat(sqrt(diag(D))',n,1);
            % Euclidean Distance
    %         v = dot(X,X,1);
    %         D = bsxfun(@plus,v,v')-2*(X'*X);
            % Hamming Distance
            D = X'*X + (1-X)'*(1-X);
            D = (size(X,1) - D)./size(X,1);
            Xbar = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                [~,idx] = min(sum(D(:,class==class_idx),2));
                Xbar(:,class_idx) = double(X(:,idx));
            end
            Xbarbar = mean(Xbar,2);
            Bi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                Bi(:,:,class_idx) = (Xbar(:,class_idx) - Xbarbar)*((Xbar(:,class_idx) - Xbarbar)');
            end
            B = sum(Bi,3);
            clear Bi
            Wi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                n = sum(class == class_idx);
                Wi(:,:,class_idx) = (1/(n-1))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))'));
            end
            W = sum(Wi,3);
            clear Wi
            [V,D] = eig(diag(diag(W))\B);
            [~,I] = max(abs(real(diag(D))));
            d = V(:,I);
            
        case 'medianLDA'
            X = double(X);
            Xbar = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                Xbar(:,class_idx) = median(X(:,class==class_idx),2);
            end
            Xbarbar = median(Xbar,2);
            Bi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                Bi(:,:,class_idx) = (Xbar(:,class_idx) - Xbarbar)*((Xbar(:,class_idx) - Xbarbar)');
            end
            B = sum(Bi,3);
            clear Bi
            Wi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                n = sum(class == class_idx);
                Wi(:,:,class_idx) = (1/(n-1))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))'));
            end
            W = sum(Wi,3);
            clear Wi
            [V,D] = eig(W\B);
            [~,I] = max(abs(real(diag(D))));
            d = V(:,I);

        case 'medianBayes'
            X = double(X);
            Xbar = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                Xbar(:,class_idx) = median(X(:,class==class_idx),2);
            end
            Xbarbar = median(Xbar,2);
            Bi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                Bi(:,:,class_idx) = (Xbar(:,class_idx) - Xbarbar)*((Xbar(:,class_idx) - Xbarbar)');
            end
            B = sum(Bi,3);
            clear Bi
            Wi = zeros(size(X,1),size(X,1),n_class);
            for class_idx = 1:n_class
                n = sum(class == class_idx);
                Wi(:,:,class_idx) = (1/(n-1))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))*...
                    ((X(:,class==class_idx) - repmat(Xbar(:,class_idx),1,n))'));
            end
            W = sum(Wi,3);
            clear Wi
            [V,D] = eig(diag(diag(W))\B);
            [~,I] = max(abs(real(diag(D))));
            d = V(:,I);

        case {'independantBernoulli' 'independantBernoulli_Bayes'}
            P = zeros(size(X,1),n_class);
            for class_idx = 1:n_class
                P(:,class_idx) = mean(X(:,class==class_idx),2);
            end
            if n_class == 2
                d = P(:,2) - P(:,1);
            else
                d = P;
            end
            
        case 'DIPPSprojection'
            if n_class ~= 2
                error('Method: "DIPPSprojection" not supported for number of clusters other than 2.')
            end
            d = DIPPSanalysis(X(:,class>0),class(class>0)==2);
            
        otherwise
            error(['Invalid method "' method '" for maldiDA()'])
            
    end
    
    if nargin == 4
        X = Xtest;
        clear Xtest
    end
    
    
    % Test
    switch method
        
        case 'independantBernoulli'
            X = double(X);
            F = zeros(n_class,size(X,2));
            for class_idx = 1:n_class
                F(class_idx,:) = sum(log(abs(repmat(P(:,class_idx),1,size(X,2)) + X - 1)));
            end
            [~,class_out] = max(F);
            
        case 'independantBernoulli_Bayes'
            X = double(X);
            if n_class ~= 2
                error('Method: "independantBernoulli_Bayes" not supported for number of clusters other than 2.')
            end
            F = zeros(n_class,size(X,2));
            for class_idx = 1:n_class
                F(class_idx,:) = sum(log(abs(repmat(P(:,class_idx),1,size(X,2)) + X - 1)));
            end
            class_out = ((F(1,:) - F(2,:)) < log(sum(class==2)/sum(class==1))) + 1;

        case 'DIPPSprojection'
            class_out = d'*X;
            class_out(class_out>0) = 2;
            class_out(class_out<0) = 1;
            
        otherwise
            if sum(isnan(d)) > 0
                class_out = nan(1,size(X,2));
            else
                if n_class == 2
                   b = d'*Xbar;
                   if b(1) > b(2)
                       d = -1*d;
                   end
                   b = mean(d'*Xbar);
                   F = d'*X;
                   class_out = (F > b)+1;
                else
                    F = zeros(n_class,size(X,2));
                    for class_idx = 1:n_class
                        F(class_idx,:) = d'*(X - repmat(Xbar(:,class_idx),1,size(X,2)));
                    end
                    [~,class_out] = min(abs(F));
                end
            end
    end
    
    valid_vars = double(valid_vars);
    valid_vars(valid_vars > 0) = d;
    d = valid_vars;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end

    