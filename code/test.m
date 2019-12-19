%% Information
% Author: Mingtao Huang
% Date  : 2019/12/19
% Brief : Test the algorithm proposed in "Permuted Sparse Representation for 3D Point Clouds"

%% Basic
clear;
close all;
% delete the followed commentary if first time running it
% mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut

%% Read .off file
% file_path = '../../../ModelNet10/ModelNet10/';
% file_name = [file_path, 'toilet/train/toilet_0001.off'];
file_name = '../data/toilet_0001.off';
[vertex,face] = read_off(file_name);

%% Generate points cloud
N = 1024;
pc = samplepc(vertex,face,N);
% scatter3(pc(1,:),pc(2,:),pc(3,:));

%% Optimization
X = pc';
D = dctmtx(N);
P = eye(N);
lambda = 0.1;
s = 0.8;
P_ = P;

for m = 1:100
    for n = 1:10
        % C-sub problem
        Y = D'*P*X;
        C = Y;
        C(abs(Y)<=lambda) = 0;
        C(Y>lambda) = C(Y>lambda) - lambda;
        C(Y<lambda) = C(Y<lambda) + lambda;

        % P-sub problem
        Xhat = D*C;
        Z = Xhat*X';
        Z = Z-min(Z(:));
        [~,P] = auction_algorithm(Z);
        fprintf('lambda: %f, diff: %d\n',lambda,full(sum(abs(P(:)-P_(:)))))
        if sum(abs(P(:)-P_(:))) < 50
            break
        end
        P_ = P;
    end
    lambda = lambda * s;
    fprintf('error: %f\n',norm(Xhat-X,'fro')/norm(X,'fro'));
    if norm(Xhat-X,'fro') < norm(X,'fro')*5e-3
        break
    end
end


