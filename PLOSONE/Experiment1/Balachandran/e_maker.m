function e = e_maker(X,Y,W) 
W1 = W(:,1,:); W2 = W(:,2,:); W3 = W(:,3,:); 
D = Y-X; 
D1 = D(1,:); D2 = D(2,:); D3 = D(3,:); 
D = reshape([repmat(D1,3,1);repmat(D2,3,1);repmat(D3,3,1)],size(W)); 
D1 = D(:,1,:); D2 = D(:,2,:); D3 = D(:,3,:); 
e = [W1.*D1 + W2.*D2 + W3.*D3]; 
e = e(:);
end