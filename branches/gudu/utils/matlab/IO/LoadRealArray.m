% Load real array  into matlab
function A = LoadRealArray(filename,IOmethod,byteorder)
% filename = '../dsigma.bin';
% CHOOSE BETWEEN byte storage format below (uncomment the one needed)
% byteorder = 'ieee-be'; % IEEE Big-Endian format 
%byteorder = 'ieee-le'; % IEEE Little-Endian format
switch IOmethod
    case 1 % unformatted file
        fid=fopen(filename,'r',byteorder);
        [nrec,count] = fread(fid,1,'int32');
        r1 = fread(fid, 1, 'int32');
        r2 = fread(fid, 1, 'int32');
        A = fread(fid, [r1 r2], 'float64');
        fclose(fid);            
    case 2 % unformatted file from ftn95
        fid=fopen(filename,'r',byteorder);
        [nrec,count] = fread(fid,1,'int8');
        [nrec,count] = fread(fid,1,'int32');
        r1 = fread(fid, 1, 'int32');
        r2 = fread(fid, 1, 'int32');
        A = fread(fid, [r1 r2], 'float64');
        fclose(fid);
    case 3 % unformatted file from gfortran (MacOS)
        fid=fopen(filename,'r',byteorder);
        [nrec,count] = fread(fid,1,'int32');
        r1 = fread(fid, 1, 'int32');
        r2 = fread(fid, 1, 'int32');
        A = fread(fid, [r1 r2], 'float64');
        fclose(fid);                
    otherwise  % binary file
        fid=fopen(filename,'r',byteorder);
        r1 = fread(fid, 1, 'int32');
        r2 = fread(fid, 1, 'int32');
        A = fread(fid, [r1 r2], 'float64');
        fclose(fid);
end
A = sparse(A);
return
