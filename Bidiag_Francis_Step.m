function [B_next, bulge] = bidiag_francis_step(B)
% Performs one step of the bidiagonal Francis algorithm to zero out the
% off-diagonal elements of the bidiagonal matrix B, using a bulge-chasing
% strategy that updates only the nonzero elements.
%
% Inputs:
%   - B: the bidiagonal matrix to be updated
%
% Outputs:
%   - B_next: the updated bidiagonal matrix
%   - bulge: the position of the new bulge (the index of the column that
%     contains the new nonzero off-diagonal element)
%
% Author: Rodrigo de Oliveira Neto
% Date: 24 April 2023

% Get the size of the matrix B
[m, n] = size(B);

% Find the first nonzero off-diagonal element
[k, l] = find(B(1:min(m-1,n),2:min(m,n)-1), 1);

if isempty(k)
    % If there are no more nonzero off-diagonal elements, we're done
    B_next = B;
    bulge = 0;
    return
end

% Compute the Givens rotation matrix that zeros out the bulge
c = B(k, k) / norm(B(k:k+1, l));
s = -B(k+1, k) / norm(B(k:k+1, l));
G = [c s; -s c];

% Apply the Givens rotation to the right of the bulge
B(k:k+1, l:n) = G' * B(k:k+1, l:n);

% Apply the Givens rotation to the left of the bulge
B(1:k+1, k:k+1) = B(1:k+1, k:k+1) * G;

% Find the new bulge
if l < n
    % If the current bulge is not at the end of the matrix, update the
    % corresponding element and set the new bulge position
    B(k, l+1) = 0;
    bulge = l+1;
else
    % Otherwise, the new bulge is in the same column as the old one
    bulge = l;
end

% Compute the Givens rotation matrix that zeros out the new bulge
c = B(k, bulge) / norm(B(k:k+1, bulge));
s = -B(k+1, bulge) / norm(B(k:k+1, bulge));
H = [c -s; s c];

% Apply the Givens rotation to the left of the new bulge
if bulge < n
    B(k:k+1, bulge+1:n) = H' * B(k:k+1, bulge+1:n);
end

% Apply the Givens rotation to the right of the new bulge
B(1:k+1, bulge:bulge+1) = B(1:k+1, bulge:bulge+1) * H';

% Return the updated matrix and the new bulge position
B_next = B;
end