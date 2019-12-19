%% Information
% Author: Mingtao Huang
% Date  : 2019/12/19
% Brief : Test the algorithm proposed in "Permuted Sparse Representation for 3D Point Clouds"

%% Basic
clear;
close all;
% delete the followed commentary if first time running it
% mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
warning('off');

%% Read .off file
% file_path = '../../../ModelNet10/ModelNet10/';
% file_name = [file_path, 'toilet/train/toilet_0001.off'];
file_name = '../data/bed_0001.off';
[vertex,face] = read_off(file_name);

%% Generate points cloud
N = 1024;
pc = samplepc(vertex,face,N);
figure();
scatter3(pc(1,:),pc(2,:),pc(3,:));

%% Optimization
X = pc';
D = dctmtx(N);
P = eye(N);
lambda = 10;
s = 0.92;
P_ = P;

error = zeros(1,100);
rate = zeros(1,100);

for m = 1:100
    for n = 1:100
        % C-sub problem
        Y = D'*P*X;
        C = Y;
        C(abs(Y)<=lambda) = 0;
        C(Y>lambda) = C(Y>lambda) - lambda;
        C(Y<-lambda) = C(Y<-lambda) + lambda;

        % P-sub problem
        Xhat = D*C;
        Z = Xhat*X';
        Z = Z-min(Z(:))+1e-5;
        [~,P] = auction_algorithm(Z);
%         fprintf('lambda: %f, diff: %d\n',lambda,full(sum(abs(P(:)-P_(:)))));
        if sum(abs(P(:)-P_(:))) < 30
            break
        else
            P_ = P;
        end
    end
    lambda = lambda * s;
    fprintf('m:%d,  n:%d\n',m,n);
    fprintf('lambda: %f, diff: %d\n',lambda,full(sum(abs(P(:)-P_(:)))));
    fprintf('compress rate: %f\n',sum(C(:)<1e-5)/numel(C));
    fprintf('error: %f\n',norm(Xhat-P*X,'fro')/norm(X,'fro'));
    error(m) = norm(Xhat-P*X,'fro')/norm(X,'fro');
    rate(m) = sum(C(:)<1e-5)/numel(C);
    if norm(Xhat-P*X,'fro') < norm(X,'fro')*3e-3
        break
    end
end

figure();
plot(1-rate(1:m), error(1:m));
xlabel('compress rate');
ylabel('error');

