ndof = 12800;
nparams = 2;
isnp = 300;
nmodes = 60;
ndeim = 60;

%---------------------------
% POD of Q

disp("POD of Q");
fflush(stdout);

X = 0.0;

filename = '../step1/param1/snapshot_u';
A=load(filename);

ii = 0;
for j = 1:isnp
for i = 1:ndof
ii = ii + 1;
X(i,j) = A(ii);
end
end

filename = '../step1/param2/snapshot_u';
A=load(filename);

ii = 0;
for j = 1:isnp
for i = 1:ndof
ii = ii + 1;
X(i,isnp+j) = A(ii);
end
end

%------------------------------------

[V,S,U] = svd(X);

for i = 1:nparams*isnp
if(S(i,i)<0.0)
disp("error in POD of Q");
fflush(stdout);
end
end

qm = V(1:ndof,1:nmodes);

% check if svd is correct
XX = V*S*(U');
zero = 0;
for j = 1:nparams*isnp
for i = 1:ndof
zero = (X(i,j)-XX(i,j))^2;
end
end
zero = sqrt(zero);
disp('zero-svd'), disp(zero);


% check for orthonormality
iden = (qm')*qm;
s1 = 0; s2 = 0;
for i = 1:nmodes
for j = 1:nmodes
if(i==j)
 s1 = s1 + iden(i,j);
else
 s2 = s2 + iden(i,j);
end
end
end

disp("one"), disp(s1/nmodes);
disp("zero"), disp(s2);

%---------------------------
% POD of F

disp("POD of F");
fflush(stdout);

X = 0.0;

filename = '../step1/param1/snapshot_rhs';
A=load(filename);

ii = 0;
for j = 1:isnp
for i = 1:ndof
ii = ii + 1;
X(i,j) = A(ii);
end
end

filename = '../step1/param2/snapshot_rhs';
A=load(filename);

ii = 0;
for j = 1:isnp
for i = 1:ndof
ii = ii + 1;
X(i,isnp+j) = A(ii);
end
end


[V,S,U] = svd(X);

for i = 1:nparams*isnp
if(S(i,i)<0.0)
disp("error in POD of F");
fflush(stdout);
end
end

UF = V(1:ndof,1:ndeim);
%--------------------------
% DEIM

disp("DEIM");
fflush(stdout);

z1 = 0;
maxval = 0.0;
for i = 1:ndof
if(abs(UF(i,1))>maxval)
maxval = abs(UF(i,1));
z1 = i;
end
end

phi = 0.0;
P = 0.0;
r = 0.0;
pvector = 0;

phi(1:ndof,1) = UF(1:ndof,1);
P(1:ndof,1) = 0.0; P(z1,1)=1.0;
pvector(1) = z1;

for l = 2:ndeim
PTphi = (P')*phi;
PTphil = (P')*UF(1:ndof,l);
alpha = inv(PTphi)*PTphil;
r = UF(1:ndof,l) - phi*alpha;

zl = 0;
maxval = 0.0;
for i = 1:ndof
if(abs(r(i))>maxval)
maxval = abs(r(i));
zl = i;
end
end

phi(1:ndof,l) = UF(1:ndof,l);
P(1:ndof,l) = 0.0; P(zl,l)=1.0;
pvector(l) = zl;

end

B_deim = (qm')*UF*inv((P')*UF);

%--------------------------
% compute qtil0

filename = '../step1/param1/init';
q0=load(filename);
qtil0 = (qm')*q0;

%--------------------------
% Save data

% qm
fileID = fopen('pod_basis','w');
for i = 1:ndof
for j = 1:nmodes
fprintf(fileID,'%24.16f \n',qm(i,j));
end
end
fclose(fileID);


% B_deim
fileID = fopen('B_DEIM','w');
for i = 1:nmodes
for j = 1:ndeim
fprintf(fileID,'%24.16f \n',B_deim(i,j));
end
end
fclose(fileID);


% pvector
fileID = fopen('pvector','w');
for i = 1:ndeim
fprintf(fileID,'%i \n',pvector(i));
end
fclose(fileID);


% P
fileID = fopen('P_DEIM','w');
for i = 1:ndof
for j = 1:ndeim
fprintf(fileID,'%24.16f \n',P(i,j));
end
end
fclose(fileID);


% qtil0
fileID = fopen('qtil0','w');
for i = 1:nmodes
fprintf(fileID,'%24.16f \n',qtil0(i));
end
fclose(fileID);




%--------------------------
