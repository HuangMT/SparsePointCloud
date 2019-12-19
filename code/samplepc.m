function pc = samplepc(vert,face,N)
% This function is to sample points from the ver-face structure to generate
% points cloud.
% Sample unifromly in each face, and the probability of each face is
% proportional to the area of face.
% -----------------------------------------
% author: Mingtao Huang
% version: 0.1.0
% last change date: 2019/12/19
% -----------------------------------------
% INPUT:
%   vert: vertex, a matrix with size 3xR, each col is the coo of a vertex.
%   face: a matrix with size 3xK, each col is 3 id, means the corresponding
%       vertex form a surface.
%   N: sample points. (default 1024)

K = size(face,2);

if nargin < 3
    N = 1024;
end

A = vert(:,face(2,:))-vert(:,face(1,:)); % the 1-st side
B = vert(:,face(3,:))-vert(:,face(1,:)); % the 2-nd side
AxB = cross(A,B,1);
area = sqrt(sum(AxB.^2,1));
prob = area./sum(area);

sfID = datasample(1:K, N, 'Weights', prob); % determin which face to choose
ab = rand(2,N);
inv_idx = sum(ab,1) > 1;
ab(:,inv_idx) = 1-ab(:,inv_idx); % determin the position in surface

pc = ab(1,:).*A(:,sfID)+ab(2,:).*A(:,sfID)+vert(:,face(1,sfID));

end