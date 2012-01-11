close all
% clear all
clc

IOmethod = 1; %0:binary ; 1:classical unformatted ; 2:unformatted ftn95

% CHOOSE BETWEEN byte storage format below (uncomment the one needed)
% byteorder = 'ieee-be'; % IEEE Big-Endian format 
byteorder = 'ieee-le'; % IEEE Little-Endian format

% Preconditioning matrix
filename = '../SparseMatrix.bin'
[Asp,coo.i,coo.j,coo.a] = LoadSparseMatrix(filename,IOmethod,byteorder);
coo.nrow = size(Asp,1);
coo.nnz  = nnz(Asp);

% Let make a test
if (0) 
    Asp = speye(3); Asp(1,2) = 4; Asp(2,1) = -5;
    [coo.i,coo.j,coo.a] = find(Asp);
    coo.j([2 4]) = coo.j([4 2]);
    coo.a([2 4]) = coo.a([4 2]);
    coo.nrow = size(Asp,1);
    coo.nnz = nnz(Asp);
    csr = ConvertSparseCOO2CSR(coo);

    b = zeros(csr.nrow,1);
    for i = 1 : csr.nrow
    disp(i)
    e = b;
    e(i) = 1;
    out = SparseCSRGaussSeidel(e,b,csr);
    out = sparse(out);
    Gsp(:,i) = out;
    end

    [M,N] = GaussSeidelsplitting(Asp,0);
    G = M\N;

end

% return

spy(Asp)

csr = ConvertSparseCOO2CSR(coo);

b = zeros(csr.nrow,1);
for i = 1 : csr.nrow
    disp(i)
    e = b;
    e(i) = 1;
out = SparseCSRGaussSeidel(e,b,csr);
    out = sparse(out);
    Gsp(:,i) = out;
end

[M,N] = GaussSeidelsplitting(Asp,0);
G = M\N;



return

% Matrix-vector product routine
filename = '../A.bin'
Adp = LoadRealArray(filename,IOmethod, byteorder);
Adp = sparse(Adp);

% load test
% Adp=A2;

subplot(2,2,1)
spy(Asp)
title('Sparse matrix - Precond Assembly')

subplot(2,2,2)
spy(Adp)
title('Sparse matrix - Direct product assembly')

% Adp = not(not(Adp));
% Asp = not(not(Asp));

B = Asp-Adp;
[i,j,s]= find(B);
idx = find(abs(s)<5e-14);
i(idx) = [];
j(idx) = [];
s(idx) = [];
B = sparse(i,j,s);
disp('elements in B discarded (magnitudes less than 5e-14)')

subplot(2,2,3)
spy(B)
title('Sparse matrix - Difference matrix')


B = full(B);

subplot(2,2,4)
mesh(B)
title('Difference matrix')

[i,j,s] = find(sparse(B));
if not(isempty(i))
    Gidx = min(i);
    disp(Gidx)

% figure
% spy(B)

Nx = 5+2;
Ny = 5+2;
Nz = 5+1;

% determine indexes automatically
k = mod(mod(Gidx,Nx*Nz),Nz);
i = mod(Gidx-k,Nx*Nz)/Nz+1;
j = (Gidx - k - (i-1)*Nz)/(Nx*Nz)+1;

% i = 2;
% j = 1;
% k = 1;

Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz;

disp(Gidx)
disp(sprintf('Coordinates for first element that are different from zero'))
disp(sprintf('in difference matrix is (i,j,k)=(%d,%d,%d)',i,j,k))

x = linspace(0,1,Nx);
y = linspace(0,1,Ny);
z = linspace(0,1,Nz);

[X,Y,Z] = meshgrid(x,y,z);
else
    disp('Preconditioning matrix matches the direct matrix-vector assembly! :o)')
end




