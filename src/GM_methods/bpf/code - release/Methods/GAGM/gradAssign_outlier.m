function [X,nbMatVec] = gradAssign(W, E12, b0, bStep, bMax, tolB, tolC, conf1, conf2)
%{
Timothee Cour, Praveen Srinivasan
implementation of graduated assignment
%}

% Xga = gradAssign(W, E12, 0.5, 1.075,10,1e-3,1e-3, target);
% Xb = bistocNormalize(X, 1e-3, 1000);

W=tril(W);
sumW=2*sum(sum(W))-sum(diag(W));
if sumW==0
    sumW=eps;
end
W=W*length(W)/sumW;

% W=trilW2W(W);

    target=[];

% discretisation=@discretisationMatching_max;
discretisation=@computeDiscreteCorrespondancesGreedy;


[n1,n2]=size(E12);
indValidMatches=find(E12);
%TODO:handle E12 not full

% % conversion
% [target,ignore]=find(target);
% target=target';

X = ones(n1,n2) + eps;
%     X0 = ones(length(indValidMatches),1) + eps;

% get groups for bistochastic normalization
[idx1 ID1] = make_groups(conf1);
[idx2 ID2] = make_groups(conf2);
if ID1(end) < ID2(end)
    [idx1 ID1 idx2 ID2 dumVal dumSize] = make_groups_slack(idx1, ID1, idx2, ID2);
    dumDim = 1;
elseif ID1(end) > ID2(end)
    [idx2 ID2 idx1 ID1 dumVal dumSize] = make_groups_slack(idx2, ID2, idx1, ID1);
    dumDim = 2;
else
    dumDim = 0; dumVal = 0; dumSize = 0;
end
idx1 = idx1-1; idx2 = idx2-1;

% maxBIters = 1000;
maxBIters = 1000;
%nthIter = round(maxIters / 10);
nthIter = 200;


isDisplay=0;
if isDisplay
    fig=figure2;
    display_result(X);
end

b=b0;
i = 1;
aIter = 1;
nbMatVec=0;
nMatch = length(W);

while ( b < bMax )

    oldX{1} = X;
    oldX{2} = X;
    %     oldM = X;
    bIter = 1;
    while (1)
        nbMatVec=nbMatVec+1;

%         X(indValidMatches) = exp(min(b * mex_w_times_x_symmetric_tril(X(indValidMatches),W), 700));
        X(indValidMatches) = exp(b * mex_w_times_x_symmetric_tril(X(indValidMatches),W));
        %         X = exp(b*( W*X(:) ));
        X=reshape(X, n1,n2);

        if any(isnan(X(:)) | isinf(X(:)) )
            X=oldX{1};
            break;
        end
        
        X2 = X(indValidMatches);
        X_slack = [X2; dumVal*ones(dumSize,1)];
        X_slack = mexBistocNormalize_match_slack(X_slack, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC, dumDim, dumVal, int32(1000));
%         pause
        X2 = X_slack(1:nMatch);   
        X(indValidMatches) = X2;
        
        if any(isnan(X(:)) | isinf(X(:)) )
            X=oldX{1};
            break;
        end


        if isDisplay
            if n1~=n2
                display_result(Xslack);
            else
                display_result(X);
            end
        end
        
%         sum(X(:))
%         pause;
        %prevent oscillations
        %         diff = sum(sum(abs(X-oldM)));
        %         oldM= X;
%         diff1 = sum(abs(X(:)-oldX{1}(:)));
%         diff2 = sum(abs(X(:)-oldX{2}(:)));
        diff1 = sum((X(:)-oldX{1}(:)).^2);
        diff2 = sum((X(:)-oldX{2}(:)).^2);
        diff=min(diff1,diff2);
        oldX{2}= oldX{1};
        oldX{1}= X;

        if (mod(i, nthIter) == 0)
            %             disp(sprintf('Diff on iter %d is %f; b is %f', i, diff, b));
            if ~isempty(target)
                Xd=feval(discretisation,X,E12);
                errorMatching=computeErrorMatching(Xd,target,E12);
                disp(sprintf('errorRate = %f', errorMatching.errorRate));
                0;
                %             [y,bestIdx] = max(X);
                %             disp(sprintf('Diff on iter %d is %f; b is %f, score is %f', i, diff, b, 1-sum([bestIdx] == target)/n1));
            end
        end

        % disp(diff);
        if ( diff < tolB  || bIter  > maxBIters)
            if diff < tolB
                % disp('diff < tolB');
            else
                % disp('bIter  > maxBIters');
            end
            break
        end
        bIter = bIter + 1;
        i = i + 1;
    end
    nbIters(aIter)=bIter;
    aIter = aIter + 1;
    b = b*bStep;

    %  bIter
end
% b
% aIter

0;



function [X,Xslack]=bistocNormalize_slack(X,tolC);
[n1,n2]=size(X);
if n1~=n2
    Xslack=X;
    if n1>n2
        Xslack(:,n2+1:n1)=1;
    else
        Xslack(n1+1:n2,:)=1;
    end
    Xslack = bistocNormalize(Xslack,tolC,1000);
    X=Xslack(1:n1,1:n2);
else
    Xslack=X;
    X = bistocNormalize(X,tolC,1000);
end



function display_result(X,isFirst);
if nargin<2
    isFirst=0;
end
%             figure(fig);
imagesc(X);
% if isFirst
%     imagesc(X);
% else
%     set(gca,'XData',X);
% end
caxis([0,1]);
pause(0.01);
disp2('bounds(X(:))');
% colormap(gray);
0;
%             pause;
