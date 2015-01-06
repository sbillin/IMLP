function C = C_maker(X,W) 
W1 = W(:,1,:); W2 = W(:,2,:); W3 = W(:,3,:); 
X1 = X(1,:); X2 = X(2,:); X3 = X(3,:); 
X = reshape([repmat(X1,3,1);repmat(X2,3,1);repmat(X3,3,1)],size(W)); 
X1 = X(:,1,:); X2 = X(:,2,:); X3 = X(:,3,:); 
C = [-W2.*X3+W3.*X2, +W1.*X3-W3.*X1, -W1.*X2+W2.*X1, W1, W2, W3]; 
C = permute(C,[1,3,2]); 
C = reshape(C,[],6);
end