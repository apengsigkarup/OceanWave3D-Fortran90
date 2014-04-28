% Load Sparse Matrix into matlab
function [A,i,j,s] = LoadSparseMatrix(filename, IOmethod, byteorder)
% filename = '../SparseMatrix.bin';
switch IOmethod
    case 1 % unformatted file
        fid=fopen(filename,'r',byteorder);
        [nrec,count] = fread(fid,1,'int32');
        nnz = fread(fid, 1, 'int32');
        i = fread(fid, nnz, 'int32');
        j = fread(fid, nnz, 'int32');
        s = fread(fid, nnz, 'float64');
        fclose(fid);              
    case 2 % unformatted file from ftn95
        fid=fopen(filename,'r',byteorder);
        [nrec,count] = fread(fid,1,'int8');
        [nrec,count] = fread(fid,1,'int32');
        nnz = fread(fid, 1, 'int32');
        i = fread(fid, nnz, 'int32');
        j = fread(fid, nnz, 'int32');
        s = fread(fid, nnz, 'float64');
        fclose(fid);
    case 3 % unformatted file from gfortran (MacOS)
        fid=fopen(filename,'r',byteorder);
        [nrec,count] = fread(fid,1,'int32');
        nnz = fread(fid, 1, 'int32');
        i = fread(fid, nnz, 'int32');
        j = fread(fid, nnz, 'int32');
        s = fread(fid, nnz, 'float64');
        fclose(fid);        
    otherwise  % binary file
        fid=fopen(filename,'r',byteorder);
        nnz = fread(fid, 1, 'int32');
        i = fread(fid, nnz, 'int32');
        j = fread(fid, nnz, 'int32');
        s = fread(fid, nnz, 'float64');
        fclose(fid);
end
A = sparse(i,j,s);
