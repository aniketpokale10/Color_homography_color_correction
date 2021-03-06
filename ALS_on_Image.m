% This code implements the Alternating least squares algorithm for color
% homography color correction.

input = imread('Input_data/macbeth_input1.jpg'); % input image is the image to be corrected
original = imread('Input_data/macbeth_original.jpg'); %%original image is how the corrected image should be 
original = rgb2xyz(original); % We take XYZ color space as reference for mapping
iterations = 50; % The alternating least squares algorithm converges in about 15 iterations 

P = []; Q = []; B = []; % convert images to array of n * 3
for i=1:size(input,1)
    for j=1:size(input,2)
        P = [P; input(i,j,1), input(i,j,2), input(i,j,3)];
        Q= [Q; original(i,j,1), original(i,j,2), original(i,j,3)];
    end
end

B = double(B);
P = double(P);
P_init = P;
n = length(P);

% ALS algorithm
D = speye(n);
for k=1:50
    for j=1:n
        if (norm(P(j,:))^2) ~= 0
            D(j,j) = (P(j,:)*Q(j,:)')/(norm(P(j,:))^2); % calculating matrix of shading factors
        end
    end

    M = pinv(D*P_init)*Q;
    P = P_init*M;        % The corrected RGB values are stored in P

end

P_x = rgb2xyz(P);
CCerrors(P_x,Q); % output the colr difference errors

P = uint8(P);
count = 1;
output = zeros(size(input)); % output image for visualization
for i=1:size(input,1)
    for j=1:size(input,2)
        output(i,j,:) = P(count,:); % converting P back to an image for visualization
        count = count+1;
    end
end

original = xyz2rgb(original,'OutputType','uint8');
subplot(2,2,3); imshow(uint8(original)); title('Original');
subplot(2,2,1);
imshow(input); title('Intput image');
subplot(2,2,2);
imshow(uint8(output)); title('Corrected image');




