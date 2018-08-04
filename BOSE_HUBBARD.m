function [HAMILTONIAN] = BOSE_HUBBARD(A_dag,A,N,A_id,J,U,mu,dim,sites,PBC)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Ham=zeros(dim,dim);

AdA=kron(A_dag,A);
AAd=kron(A,A_dag);

for i=(0:sites-2)
    
    ID1=A_id{i+1};
    ID2=A_id{sites-1-i};
    Ham=Ham+J*(kron(ID1,kron(AdA,ID2))+kron(ID1,kron(AAd,ID2)));
end

N2=N*N; 
for i=(0:sites-1)
   ID1=A_id{i+1};
   ID2=A_id{sites-i};
   Ham=Ham+0.5*U*(kron(ID1,kron(N2,ID2)))+mu*(kron(ID1,kron(N,ID2))) ;
end

Ham=Ham+PBC*J(kron(A_dag,kron(A_id{sites-2},A))+kron(A,kron(A_id{sites-2},A_dag)));
HAMILTONIAN=Ham;