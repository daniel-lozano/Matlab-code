function [HAMILTONIAN] = Hamiltonian_t(N,dim,Dx,Dy,Dz,hx,hy,hz,A,t,HOMOG,STAG,RANDOM,PBC)
% Halla el hamiltoniano en cuestion dependiendo de los parametros
% ingresados

%matrices de Pauli

sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
xx=kron(sx,sx); 
yy=kron(sy,sy); 
zz=kron(sz,sz);

Ham=zeros(dim,dim);

%Caso periodico
Case_pbc=Dx*kron(sx,kron(eye(2.^(N-2)),sx)) ...;
    +Dy*kron(sy,kron(eye(2.^(N-2)),sy)) ...;
    +t*Dz*kron(sz,kron(eye(2.^(N-2)),sz));

%Terminos del Hamiltoniano sin campo 
for i=(0:N-2)
    
    Id1=A{i+1};
    Id2=A{N-1-i};  
    Ham=Ham+Dx*kron(Id1,kron(xx,Id2)) ... 
        +Dy*kron(Id1,kron(yy,Id2)) ... 
        +t*Dz*kron(Id1,kron(zz,Id2))   ;
    
end

%Terminos del Hamiltoniano con campo 

%Caso de campo homogeneo
if (HOMOG==1 && STAG==0 && RANDOM==0)
    
for i=(0:N-1)
    Id1=A{i+1}; 
    Id3=A{N-i};
    Ham=Ham+hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3)) ;
end


%Caso de campo alternado
elseif (HOMOG==0 && STAG==1 && RANDOM==0)
    
for i=(0:N-1)
    Id1=A{i+1};  
    Id3=A{N-i};
    Ham=Ham+ ((-1.0).^(i)) *(hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3))) ;
end

%Caso de campo aleatorio
elseif (HOMOG==0 && STAG==0 && RANDOM==1)
    
    
    rng('shuffle');%cambia la semilla a una dependiente del tiempo
    
    
for i=(0:N-1)
    Id1=A{i+1};
    Id3=A{N-i};
    Ham=Ham+ (2*rand()-1) *(hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3))) ;
end

else
    disp("No es claro que caso de campo es necesario")
end






HAMILTONIAN = Ham+PBC*Case_pbc;

end

