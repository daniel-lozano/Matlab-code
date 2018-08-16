function [HAMILTONIAN,NUMERO] = BOSE_HUBBARD_FUNCTION(A_dag,A,N,A_id,J,U,mu,dim,dim_site,sites,PBC)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Ham=zeros(dim,dim);
Numero=zeros(dim,dim);

AdA=kron(A_dag,A);% dim 2
AAd=kron(A,A_dag);% dim 2

for i=(1:sites-1)
    ID1=A_id{i};%matriz identidad a la i-1
    ID2=A_id{sites-i};% matriz identidad a la sites-i-1
    Ham=Ham-J*(kron(ID1,kron(AdA,ID2))+kron(ID1,kron(AAd,ID2)));%dim= i-1+sites-i-1+2=sites
end


N2=N*(N-eye(dim_site)); 
for i=(1:sites)
   ID1=A_id{i};
   ID2=A_id{sites-i+1};  
   
%    disp("sitio "+num2str(i))
%    disp(size(ID1))
%    disp(size(ID2))
   
   Ham=Ham+0.5*U*(kron(ID1,kron(N2,ID2)))+mu*(kron(ID1,kron(N,ID2))) ;
   Numero=Numero+kron(ID1,kron(N,ID2));
   %disp(i-1)
end


Ham=Ham+PBC*J*(kron(A_dag,kron(A_id{sites-1},A))+kron(A,kron(A_id{sites-1},A_dag)));
HAMILTONIAN=Ham;
NUMERO=Numero;
