function y = orthonormalization(E, X)
n = length(E);
y = zeros(n);
for i = 1:n
    v = E(:,i);
    for j = 1:i-1
        proj_ij = (E(:,i)'*X*y(:,j))/(y(:,j)'*X*y(:,j));
        v = v-proj_ij* y(:,j);
    end
    y(:,i) = v/(v'*X*v)^0.5;
end