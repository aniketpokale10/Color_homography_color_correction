function H = ransac(p,q)

p = p';  q = q'; 
mit = 5000; % max number of iterations

if ~all(size(p)==size(q))
    error('Data sets must have the same dimension');
end

npts = length(p);
s = 4; % minimum number of points to solve for homography


% RANSAC loop---------------------------------------------
prob = 0.99;
maxtrials = 1000;
bestM = [];      % Sentinel value allowing detection of solution failure.
trialcount = 0;
bestscore =  0;
N = 1;            % Dummy initialisation for number of trials.
maxDataTrials = 100;
feedback = 0;
maxTrials = 1000; % maximum number of iterations for ransac

while N > trialcount
        
        % Select at random s datapoints to form a trial model, M.
        % In selecting these points we have to check that they are not in
        % a degenerate configuration.
        degenerate = 1;
        count = 1;
        while degenerate
            % Generate s random indicies in the range 1..npts
            ind = randsample(npts, s);
            
            % Test that these points are not a degenerate configuration.
            degenerate = isdegenerate (p(:,ind),q(:,ind));
            
            if ~degenerate
                % Fit model to this random selection of data points.
                % Note that M may represent a set of models that fit the data in
                % this case M will be a cell array of models
                M = uea_H_from_x_als(p(:,ind),q(:,ind));
                
                % Depending on your problem it might be that the only way you
                % can determine whether a data set is degenerate or not is to
                % try to fit a model and see if it succeeds.  If it fails we
                % reset degenerate to true.
                if isempty(M)
                    degenerate = 1;
                end
            end

            % Safeguard against being stuck in this loop forever
            count = count + 1;
            if count > maxDataTrials
                if feedback
                    warning('Unable to select a nondegenerate data set');
                end
                break
            end
        end

        if ~degenerate
            % Once we are out here we should have some kind of model...
            % Evaluate distances between points and model returning the indices
            % of elements in x that are inliers.  Additionally, if M is a cell
            % array of possible models 'distfn' will return the model that has
            % the most inliers.  After this call M will be a non-cell object
            % representing only one model.
            [inliers, M] = homogdist( M, p, q, 0.2);
        else
            inliers = [];
        end
        
        % Find the number of inliers to this model.
        ninliers = length(inliers);
        
        if ninliers > bestscore    % Largest set of inliers so far...
            bestscore = ninliers;  % Record data for this model
            bestinliers = inliers;
            bestM = M;
            
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  ninliers/npts;
            pNoOutliers = 1 - fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            N = log(1-prob)/log(pNoOutliers);
        end
        
        trialcount = trialcount+1;
        if feedback
            fprintf('trial %d out of %d         \r',trialcount, ceil(N));
        end

        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            warning('ransac reached the maximum number of %d trials',...
                    maxTrials);
            break
        end
end

H = bestM'; % final homography matrix

end



function [inliers, H] = homogdist(H, lx1, lx2, t)

    % Calculate, in both directions, the transfered points    
    Hx1 = H*lx1;
    
    white = whitepoint;
    
    % Calculate lab distance
    luv_ref = HGxyz2luv(lx2',white)'; % reference LUV
    luv_est = HGxyz2luv(Hx1',white)'; % reference LUV

    uv_ref = bsxfun(@rdivide,luv_ref(2:3,:),max(luv_ref(1,:),eps));
    uv_est = bsxfun(@rdivide,luv_est(2:3,:),max(luv_est(1,:),eps));

    d = sqrt(sum((uv_ref-uv_est).^2,1));
    inliers = find(d<t);
end



function r = isdegenerate(lx1,lx2)
   
    % generate points combinations
    xcomb = combnk(1:4,3);
    ncomb = size(xcomb,1);
    
    ir1 = arrayfun(@(i) iscolinear_n(lx1(:,xcomb(i,:))), 1:ncomb);
    ir2 = arrayfun(@(i) iscolinear_n(lx2(:,xcomb(i,:))), 1:ncomb);

    r = any([ir1,ir2]);
end

function r = iscolinear_n(P)

    if ~(size(P,1)>=3)                              
        error('points must have the same dimension of at least 3');
    end
    
	r =  norm(cross(P(:,2)-P(:,1), P(:,3)-P(:,1))) < eps;
end

function [H,err,pD] = uea_H_from_x_als(p1,p2,max_iter,tol,k)

if nargin<3, max_iter = 50; end
if nargin<4, tol = 1e-20; end
if nargin<5, k = 'lin'; end

[Nch,Npx] = size(p1);

ind1 = sum(p1>0 & p1<Inf,1)==Nch;
ind2 = sum(p2>0 & p2<Inf,1)==Nch;
vind = ind1 & ind2;

if (size(p1) ~= size(p2))
 error ('Input point sets are different sizes!')
end

% definition for Graham
P = p1;
Q = p2;
N = P;
D = speye(Npx);

errs = Inf(max_iter+1,1); % error history

% solve the homography using ALS
n_it = 1; d_err = Inf;
while ( n_it-1<max_iter && d_err>tol )
    n_it = n_it+1; % increase number of iteration

    D = SolveD1(N,Q);

    P_d = P*D;
    if strcmp(k,'lin')
        M = Q(:,vind)/P_d(:,vind); % update M
    else
        K = mean(diag(P_d*P_d'))./1e3;
        M = ((P_d*P_d'+K*eye(Nch))\P_d*(Q'))';
    end
    N = M*P;

    NDiff = (N*D-Q).^2; % difference
    errs(n_it) = mean(mean(NDiff(:,vind))); % mean square error
    d_err = errs(n_it-1) - errs(n_it); % 1 order error
end

H = M;
err = errs(n_it);

pD = D;

% plot(errs); hold on;
% fprintf('ALS %d: %f\n',n_it,errs(n_it));

% figure; imagesc(reshape(diag(D),4,6));

    function D = SolveD1(p,q)
        [nCh,nPx] = size(p);
        d = (ones(1,nCh)*(p.*q))./(ones(1,nCh)*(p.*p));
        D = spdiags(d',0,nPx,nPx);
    end

end

% HGxyz2luv converts XYZ to LUV
% Arguments:
%          xyz   - Nx3 matrix for XYZ
%          white - white point XYZ values, 1x3 vector

% Returns:
%          luv - Nx3 matrix for LUV
%          up  - u'.
%          vp  - v'.

function [luv,up,vp] = HGxyz2luv(xyz,white)

if (size(xyz,2)~=3)
   disp('xyz must be n by 3'); return;   
end
luv = zeros(size(xyz,1),size(xyz,2));

% compute u' v' for sample
up = 4*xyz(:,1)./(xyz(:,1) + 15*xyz(:,2) + 3*xyz(:,3));
vp = 9*xyz(:,2)./(xyz(:,1) + 15*xyz(:,2) + 3*xyz(:,3));
% compute u' v' for white
upw = 4*white(1)/(white(1) + 15*white(2) + 3*white(3));
vpw = 9*white(2)/(white(1) + 15*white(2) + 3*white(3));

index = (xyz(:,2)/white(2) > 0.008856);
luv(:,1) = luv(:,1) + index.*(116*(xyz(:,2)/white(2)).^(1/3) - 16);  
luv(:,1) = luv(:,1) + (1-index).*(903.3*(xyz(:,2)/white(2)));  

luv(:,2) = 13*luv(:,1).*(up - upw);
luv(:,3) = 13*luv(:,1).*(vp - vpw);

end


