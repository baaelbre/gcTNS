function y = orthonormalization(E, X)
n = length(E);
E_ = zeros(n);
for i = 1:n
    v = E(:,i);
    for j = 1:i-1
        proj(i,j) = (E(:,i)'*X*E_(:,j))/(E_(:,j)'*X*E_(:,j));
        v = v-proj(i,j)* E_(:,j);
    end
    E_(:,i) = v/(v'*X*v)^1/2;
y=E_;
end