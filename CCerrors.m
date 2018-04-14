function CCerrors(x1,x2)

% Calculating luv error
luv_est = xyz2luv(x1,whitepoint);
luv_ref = xyz2luv(x2,whitepoint);

uv_ref = bsxfun(@rdivide,luv_ref(:,2:3),max(luv_ref(:,1),eps));
uv_est = bsxfun(@rdivide,luv_est(:,2:3),max(luv_est(:,1),eps));

fprintf('Luv errors = \n');
deltaE = sqrt(sum((uv_est-uv_ref).^2, 2));

mean_error = mean(deltaE);
median_error = median(deltaE);
max_error = max(deltaE);
quantile95 = max(quantile(deltaE, 0.95));

fprintf('Mean Error: %f\n',mean_error);
fprintf('Median Error: %f\n',median_error);
fprintf('95 percernt quantile Error: %f\n',quantile95);
fprintf('Max Error: %f\n\n',max_error);


% Calculating Lab error
lab_est = xyz2lab(x1);
lab_ref = xyz2lab(x2);

ab_ref = bsxfun(@rdivide,lab_ref(:,2:3),max(lab_ref(:,1),eps));
ab_est = bsxfun(@rdivide,lab_est(:,2:3),max(lab_est(:,1),eps));

fprintf('lab errors = \n');
deltaE = sqrt(sum((ab_est-ab_ref).^2, 2));

mean_error = mean(deltaE);
median_error = median(deltaE);
max_error = max(deltaE);
quantile95 = max(quantile(deltaE, 0.95));

fprintf('Mean Error: %f\n',mean_error);
fprintf('Median Error: %f\n',median_error);
fprintf('95 percernt quantile Error: %f\n',quantile95);
fprintf('Max Error: %f\n\n',max_error);

