input = imread('macbeth_input1.jpg'); % input image is the image to be corrected
original = imread('macbeth_original.jpg'); %%original image is how the corrected image should be 

input = double(input); original = double(original);

P = []; Q = []; % convert images to array of n * 3
for i=1:size(input,1)
    for j=1:size(input,2)
        P = [P; input(i,j,1), input(i,j,2), input(i,j,3)];
        Q= [Q; original(i,j,1), original(i,j,2), original(i,j,3)];
    end
end

% Preconditioning data for ransac
mean_P = mean(P); 
Centered_P = P - repmat(mean_P,[size(P,1),1]);
var_P = var(Centered_P);
sd_P = sqrt(var_P);
Normalized_P = Centered_P./repmat(sd_P,[size(P,1),1]);

mean_Q = mean(Q);
Centered_Q = Q - repmat(mean_Q,[size(Q,1),1]);
var_Q = var(Centered_Q);
sd_Q = sqrt(var_Q);
Normalized_Q = Centered_Q./repmat(sd_Q,[size(Q,1),1]);

H = ransac(Normalized_P,Normalized_Q);

if isempty(H) 
    error('Empty homography matrix, please try again');
end

% Calculating luv error
luv_est = xyz2luv(P*H,whitepoint);
luv_ref = xyz2luv(Q,whitepoint);

% d = sqrt(sum((luv_est - luv_ref).^2,1));
% d = d/length(P) ;
% d = sum(d)/3;
% fprintf('Luv Error: %f\n',d);

mean_error = 0; % Calculate the mean error
deltaE = 0;
error_vector = [];
for i=1:length(P)
    deltaEtemp = 0;
    for k=1:3
        mean_error = mean_error + abs(luv_ref(i,k) - luv_est(i,k))./(3*length(P));
        deltaEtemp = deltaEtemp + (luv_ref(i,k) - luv_est(i,k)).^2;
        error_vector = [error_vector; abs(luv_ref(i,k) - luv_est(i,k))];
    end
    deltaE = deltaE + sqrt(deltaEtemp);
end
deltaE = deltaE/(3*length(P));
fprintf('Mean Luv Error: %f\n',mean_error);
fprintf('Median Luv Error: %f\n',median(error_vector));
fprintf('delta E Luv Error: %f\n',deltaE);

O = uint8(P*H); 
count = 1;
output = zeros(size(input));
for i=1:size(input,1)
    for j=1:size(input,2)
        output(i,j,:) = O(count,:); % converting O back to an image for visualization
        count = count+1;
    end
end

% figure;
subplot(2,2,3); imshow(uint8(original)); title('Original');
subplot(2,2,1);
imshow(uint8(input)); title('Intput image');
subplot(2,2,2);
imshow(uint8(output)); title('Corrected image');



