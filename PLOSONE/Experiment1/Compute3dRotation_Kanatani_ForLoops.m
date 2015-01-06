function [data] = Compute3dRotation_Kanatani_ForLoops(xobs,yobs,Mx,My,maxIter, threshAng)
%
% For loop implementation
%

numPts = size(xobs,1);

R = eye(3);
q = zeros(4,1);
N = zeros(4,4);

% Step 1: Init Xa
xobs_t = xobs';
yobs_t = yobs';
Xa = zeros(3,4,numPts);
for i=1:numPts
  Xa(:,1,i) = yobs_t(:,i) - xobs_t(:,i);
  Xa(:,2:4,i) = skew(yobs_t(:,i) + xobs_t(:,i));
end

% Step 2: Init c, Wa
c = 0;
Wa = repmat(eye(3),1,1,numPts);

numIter = 0;
lambda = 1;
dTheta = threshAng + 1;

vThetas = zeros(maxIter,1);
vAxis = zeros(maxIter,3);

%while ( abs(lambda) > 0.1 && numIter < maxIter )   % Kanatani's termination criteria
while ( dTheta > threshAng && numIter < maxIter)

  numIter = numIter + 1;
  
  % Step 3: M
  M = zeros(4,4);
  for i=1:numPts
    M = M + Xa(:,:,i)'*Wa(:,:,i)*Xa(:,:,i);
  end

  % Step 4: N
  Mx_p_My = Mx + My;
  Mx_m_My = Mx - My;
  n0 = 0;
  for i=1:numPts
    n0 = n0 + sum(sum(Wa(:,:,i).*Mx_p_My(:,:,i)));
  end
  n = 0;
  for i=1:numPts
    A = Wa(:,:,i) * Mx_m_My(:,:,i);
    A = (A - A')/2;
    n = n + [A(3,2) A(1,3) A(2,1)]';
  end
  n = -2*n;  
  % % check if exterior matrix product should have been here
  % %  Ans: this did not help
  % n = 0;
  % for p=1:numPts
  %   A = zeros(3,3);
  %   for i=1:3
  %     for j=1:3
  %       for k=1:3
  %         for l=1:3
  %           for m=1:3
  %             for n=1:3
  %               eps1 = EddingtonEps(i,k,l);
  %               eps2 = EddingtonEps(j,m,n);
  %               A(i,j) = A(i,j) + eps1*eps2*Wa(k,m,p)*Mx_m_My(l,n,p);
  %             end
  %           end
  %         end
  %       end
  %     end
  %   end
  %   A = (A - A')/2;
  %   n = n + [A(3,2) A(1,3) A(2,1)]';    
  % end
  % n = -2*n;      
  Np = zeros(3,3);  % N' = Sum_i(Wai x (Mxi+Myi))
  for p=1:numPts
    for i=1:3
      for j=1:3
        for k=1:3
          for l=1:3
            for m=1:3
              for n=1:3
                eps1 = EddingtonEps(i,k,l);
                eps2 = EddingtonEps(j,m,n);
                Np(i,j) = Np(i,j) + eps1*eps2*Wa(k,m,p)*Mx_p_My(l,n,p);
              end
            end
          end
        end
      end
    end
  end
  N(1,1) = n0;
  N(2:4,1) = n;
  N(1,2:4) = n';
  N(2:4,2:4) = Np;

  % Step 5: eigenvalue of Mh
  Mh = M - c*N;
  [V,D] = eig(Mh);
  lambda = min(diag(D));
  idx = find(diag(D) == lambda);
  q = V(:,idx);
  q0 = q(1);
  ql = q(2:4);
  
  % Step 6: update c & Wa
  %  Note: most of the run-time was in this step
  %        step 4 comes close 2nd
  c = c + lambda/(q'*N*q); 
  for p=1:numPts    
    Tmp1 = q0*q0*Mx_p_My(:,:,p);
    S = zeros(3,3);
    S(:,1) = cross(ql,Mx_m_My(:,1,p));
    S(:,2) = cross(ql,Mx_m_My(:,2,p));
    S(:,3) = cross(ql,Mx_m_My(:,3,p));
    Tmp2 = q0*(S+S');
    
    ql_Mx_p_My_ql = zeros(3,3);
    for i=1:3
      for j=1:3
        for k=1:3
          for l=1:3
            for m=1:3
              for n=1:3
                eps1 = EddingtonEps(i,k,l);
                eps2 = EddingtonEps(j,m,n);
                ql_Mx_p_My_ql(i,j) = ql_Mx_p_My_ql(i,j) + eps1*eps2*ql(k)*ql(m)*Mx_p_My(l,n,p);
              end
            end
          end
        end
      end
    end  
    
    Wa(:,:,p) = inv(Tmp1 - Tmp2 + ql_Mx_p_My_ql);
  end

  Rprev = R;
  R = quat2rot(q);  
  dR = R*Rprev';
  [dAxis, dTheta] = rot2AxisAngle(dR);
  dTheta = dTheta*180/pi;
  
  vTheta(numIter) = dTheta;
  vAxis(numIter,:) = dAxis';
end

data{1} = R;
data{2} = numIter;
data{3} = vTheta(1:numIter)';
data{4} = vAxis(1:numIter,:);
end

function [v] = EddingtonEps(i,j,k)
  v = (j-i)*(k-j)*(k-i)/2;
end
