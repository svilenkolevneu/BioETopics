%% Problem 2 
% partial least squares done using plsregres in matlab 
% and total least squares done using SVD (algorithm taken from wikipedia)
clear; clc;
nums = [50, 100, 500, 1000];
sigs = [0.1, 0.25, 0.5, 0.75, 1];
for i = 1:length(nums);
    for j = 1:length(sigs);
        run = totvspar(nums(i),sigs(j));
        subplot(5, 4,sub2ind([4, 5], i, j));
        hist(run,30);
        title('N='+string(nums(i))+' var='+string(sigs(j)));
        xlim([0,5]);
    end
end

function slopes = totvspar(n,sig)
slopes = zeros(1000, 2);
for i = 1:1000;
    x = rand(n,1);
    y = 2*x;
    xnos = normrnd (0, sqrt(sig), n, 1);
    ynos = normrnd (0, sqrt(sig), n, 1);
    xplus = x + xnos;
    yplus = y + ynos;
    
    [XL, YL,~,~,b] = plsregress(xplus,yplus,1);
    slopes(i,1) = b(2);
    slopes(i,2) = tls(xplus,yplus);
end
end

% total least square algorithm using SVD from wikipedia
function B = tls(X,Y)
[m n]   = size(X);            % n is the width of X (X is m by n)
Z       = [X Y];              % Z is X augmented with Y.
[U S V] = svd(Z,0);           % find the SVD of Z.
VXY     = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY     = V(1+n:end,1+n:end); % Take the bottom-right block of V.
B       = -VXY/VYY;
end