function CCerrors(x1,x2)

% Calculating luv error
luv_est = xyz2luv(x1,whitepoint);
luv_ref = xyz2luv(x2,whitepoint);

uv_ref = bsxfun(@rdivide,luv_ref(:,2:3),max(luv_ref(:,1),eps));
uv_est = bsxfun(@rdivide,luv_est(:,2:3),max(luv_est(:,1),eps));

fprintf('Luv errors = \n');

mean_error = mean(mean(abs(uv_ref - uv_est)));
median_error = median(sum(abs(uv_ref - uv_est),2)./2);
max_error = max(sum(abs(uv_ref - uv_est),2)./2);
deltaE = sum(sqrt(sum((uv_est-uv_ref).^2, 2)))/length(uv_est);
quantile95 = max(quantile(uv_ref - uv_est, 0.95));

fprintf('Mean Error: %f\n',mean_error);
fprintf('Median Error: %f\n',median_error);
fprintf('Max Error: %f\n',max_error);
fprintf('delta E Error: %f\n',deltaE);
fprintf('95 percernt quantile Error: %f\n\n',quantile95);


% Calculating Lab error
lab_est = xyz2lab(x1);
lab_ref = xyz2lab(x2);

ab_ref = bsxfun(@rdivide,lab_ref(:,2:3),max(lab_ref(:,1),eps));
ab_est = bsxfun(@rdivide,lab_est(:,2:3),max(lab_est(:,1),eps));

fprintf('lab errors = \n');

mean_error = mean(mean(abs(ab_ref - ab_est)));
median_error = median(sum(abs(ab_ref - ab_est),2)./2);
max_error = max(sum(abs(ab_ref - ab_est),2)./2);
deltaE = sum(sqrt(sum((ab_est-ab_ref).^2, 2)))/length(ab_est);
quantile95 = max(quantile(ab_ref - ab_est, 0.95));

fprintf('Mean Error: %f\n',mean_error);
fprintf('Median Error: %f\n',median_error);
fprintf('Max Error: %f\n',max_error);
fprintf('delta E Error: %f\n',deltaE);
fprintf('95 percernt quantile Error: %f\n',quantile95);

