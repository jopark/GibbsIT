function MDV = cauchy (a, b)
%3.0 authors: Jun Park and Sara Rubin, Princeton University
%1.0-2.1 author: Lake-Ee Quek, AIBN
%function MDV = cauchy (a, b);
%cauchy product to generate MID vector
%a, b: input MID vectors
%MDV: MID vector from resulting cauchy product
la=length(a);
lb=length(b);
dim = la + lb - 1;
% to stagger vector a into a matrix
%resulting matrix has dim number of rows
A_write=zeros(dim+la-1,dim);
B_write=zeros(dim,1);
for i = 1:dim;
    A_write(i:i+la-1,i) = a;
end
A_write(dim+1:end,:)=[];
B_write(1:lb,1)=b;
%resulting vector is the first column of the MDV matrix
MDV=(A_write*B_write)';