clear;
clc;

fprintf('\n####################################################################\n')
fprintf(' Findpeaks: Computational Crystallography ')
fprintf('\n####################################################################\n')
fprintf(' Kayode Olumoyin ')
fprintf(' COMS 7100 ')
fprintf(' Project 2 ')
fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long e % using the highest matlab precision

warning('off','all') % This suppresses the RCOND warnings
prompt = 'Enter the input text file: ';

filename = input(prompt,'s'); %proj2_data.txt

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(filename,'r');

header = textscan(fileID,' %s %s %s ',1,'Delimiter',',');
%celldisp(header(1,3))

fprintf('\n comment = %s', string(header(1,1)))
fprintf(', %s', string(header(1,2)))
fprintf(', %s\n', string(header(1,3)))

header2 = textscan(fileID,' %s %f %f %f %f %f %f',1);
fprintf('\n')
num_unit_cell_param = 6;

unit_cell_vec = zeros(num_unit_cell_param,1);

for k=1:num_unit_cell_param
    unit_cell_vec(k) = str2double(string(header2(1,k+1)));
    fprintf('Unit cell parameter %i = %f\n', k,str2double(string(header2(1,k+1))))
end

fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab function to detect the input file as a table 
opts = detectImportOptions(filename); 

% set table as a table of double type
file_double = readtable(filename,opts);

% check the variables detected and the datatype detected
% disp([opts.VariableNames' opts.VariableTypes'])

% set table as a table of string type
opts = setvartype(opts,'string');
file_string = readtable(filename,opts);

%preview(filename,opts)

fprintf('\n####################################################################\n')
% detect the h data
h = file_double(1:end,1);
H = table2array(h(:,1));

% detect the k data
k = file_double(1:end,2);
K = table2array(k(:,1));

% detect the l data
l = file_double(1:end,3);
L = table2array(l(:,1));

% detect the F data
f = file_double(1:end,4);
F = table2array(f(:,1));

% detect the A(hkl) data
A = file_double(1:end,5);
A_hkl = table2array(A(:,1));

% detect the B(hkl) data
B = file_double(1:end,6);
B_hkl = table2array(B(:,1));

fprintf('\n####################################################################\n')

data_len = length(H);
fprintf('\n Number of hkl data = %i \n', data_len)

fprintf('\n####################################################################\n')
a = unit_cell_vec(1);
b = unit_cell_vec(2);
c = unit_cell_vec(3);

alpha = unit_cell_vec(4);
beta = unit_cell_vec(5);
gamma = unit_cell_vec(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = zeros(3);

G(1,1) = mtimes(a,a);
G(1,2) = mtimes(a,b)*cosd(gamma);
G(1,3) = mtimes(a,c)*cosd(beta);
G(2,1) = mtimes(b,a)*cosd(gamma);
G(2,2) = mtimes(b,b);
G(2,3) = mtimes(b,c)*cosd(alpha);
G(3,1) = mtimes(c,a)*cosd(beta);
G(3,2) = mtimes(c,b)*cosd(alpha);
G(3,3) = mtimes(c,c);

fprintf('\n Metric tensor G\n')
disp(G);

Vc = sqrt(det(G));
fprintf('\n Unit cell volume, Vc computed from G = %f\n',Vc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = zeros(3);

M(1,1) = a;
M(1,2) = b*cosd(gamma);
M(1,3) = c*cosd(beta);
M(2,1) = 0;
M(2,2) = b*sind(gamma);
M(2,3) = c*(cosd(alpha) - cosd(beta)*cosd(gamma))/sind(gamma);
M(3,1) = 0;
M(3,2) = 0;
M(3,3) = c*(sqrt(sind(gamma)^2 - ...
    (cosd(alpha)^2 + cosd(beta)^2 - 2*cosd(alpha)*cosd(beta)*cosd(gamma)))/sind(gamma));

fprintf('\n Transformation matrix CRYSTAL to CARTESIAN, M matrix\n')
disp(M);

fprintf('\n Transformation matrix CARTESIAN to CRYSTAL, M_inv matrix G\n')
M_inv = inv(M);
disp(M_inv);

fprintf('\n Metric tensor G, using G = M.T*M \n')
disp(M'*M);

fprintf('\n####################################################################\n')

fprintf('\n peak search limits in FRACTIONAL COORDINATES:\n')
fprintf('0.1 <= a <= 0.9\n')
fprintf('0.1 <= b <= 0.9\n')
fprintf('0.1 <= c <= 0.9\n')

% fractional coordinates
fprintf('\n Fractional coordinates, xf \n')
xf = [0.8, 0.6, 0.4]';
disp(xf);

% Cartesian coordinates
fprintf('\n Cartesian coordinates, xc \n')
xc = mtimes(M,xf);
disp(xc);

% compute rho using fractional coordinates
rhoff = rhoxf(xf, A_hkl, B_hkl, H, K, L, Vc, data_len);
fprintf('\n Rho(xf) in fractional coordinates\n')
disp(rhoff);

% compute rho using cartesian coordinates
rhocc = rhoxc(xc, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
fprintf('\n Rho(xc) in cartesian coordinates\n')
disp(rhocc);

% compute gradient in Cartesian coordinates
grad_c = gradient_c(xc, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
fprintf('\n gradient, grad rho(xc) in cartesian coordinates\n')
disp(grad_c);

% compute gradient norm in Cartesian coordinates
fprintf('\n gradient norm, |grad rho(xc)| in cartesian coordinates = %e\n', norm(grad_c))

% Hessian in Cartesian coordinates
Hess_c = hessian_c(xc, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
fprintf('\n Hessian, hess rho(xc) in cartesian coordinates\n')
disp(Hess_c);

% eigenvalues and eigenvector of the Hessian in Cartesian coordinates
[eigVe, eigVa] = eigg(xc, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);

fprintf('\n eigenvalues of the Hessian, hess rho(xc) in cartesian coordinates\n')
disp(diag(eigVa));

fprintf('\n eigenvectors of the Hessian, hess rho(xc) in cartesian coordinates\n')
% eigVe = eigVe';
disp(eigVe);

fprintf('\n -inv(hess rho(xc))*(grad rho(xc)), h = \n')
h = linsolve(-Hess_c,grad_c);
disp(h);

fprintf('\n Cartesian coordinates xc + h \n')
new_xc = xc + h;
disp(new_xc);

new_rhocc1 = rhoxc(new_xc, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
fprintf('\n Rho(xc+h) in cartesian coordinates\n')
disp(new_rhocc1);

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%% Finding Peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xf_start = [0.1, 0.1, 0.1]';    % start search here, fractional coordinates
xc_start = mtimes(M,xf_start);  % cartesian coordinates
% disp(xc_start);

xf_end = [0.9, 0.9, 0.9]';      % end search here, fractional coordinates
xc_end = mtimes(M,xf_end);      % cartesian coordinates
% disp(xc_end);

rho_xc_start = rhoxc(xc_start, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
% disp(rhoxc0);

grad_xc_start = gradient_c(xc_start, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
Hess_xc_start = hessian_c(xc_start, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);

hh = linsolve(-Hess_xc_start,grad_xc_start);
xc_next = xc_start + hh;
%disp(xc_next);

rho_xc_next = rhoxc(xc_next, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len);
%disp(rho_xc_next);

% length in the fractional coordinates

len_a = a*(0.9 - 0.1); % length in the X direction
len_b = b*(0.9 - 0.1); % length in the Y direction
len_c = c*(0.9 - 0.1); % length in the Z direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n search using Newton-Raphson \n')
fprintf('\n length in X direction = %f angstroms\n',len_a)
fprintf('\n length in Y direction = %f angstroms\n',len_b)
fprintf('\n length in Z direction = %f angstroms\n',len_c)

% Grid spacing 0.4 angstrom

% Grid step size in fractional coordinates

GS_a = 0.4*(0.9 - 0.1)/len_a;
GS_b = 0.4*(0.9 - 0.1)/len_b;
GS_c = 0.4*(0.9 - 0.1)/len_c;

% Grid step size (FRA) along CRYSTAL coordinate axes
GS_abc_f = [GS_a, GS_b, GS_c]';
fprintf('\n Grid step size (FRA) along CRYSTAL coordinate axes = \n ')
disp(GS_abc_f');

% Grid step size in cartesian coordinates

GS_abc_c = mtimes(M,GS_abc_f); 
fprintf('\n Grid step size (CAR) along CARTESIAN coordinate axes = \n ')
disp(GS_abc_c');

% All the 3D grid points in Fractional coordinates

F_x = 0.1:GS_a:0.9+GS_a;
F_y = 0.1:GS_b:0.9+GS_b;
F_z = 0.1:GS_c:0.9+GS_c;

len_Fx = length(F_x(:));
len_Fy = length(F_y(:));
len_Fz = length(F_z(:));

fprintf('\n Number of grid points along CRYSTAL coordinate axes = \n')
grid = [len_Fx, len_Fy, len_Fz];
disp(grid);

T = len_Fx*len_Fy*len_Fz; % total length of grid points
fprintf('\n Total number of grid points\n')
disp(T);

F_x = F_x(:);
F_y = F_y(:);
F_z = F_z(:);

% search list in Fractional coordinates
llists = [];
for i = 1:len_Fx
    sublist =[];
    for j = 1:len_Fy
        ssublist=[];
        for k = 1:len_Fz
            ssublist = [ssublist; [F_x(i),F_y(j),F_z(k)]];
        end
        sublist = [sublist; ssublist];
    end
    llists = [llists; sublist];
end
list_f = reshape(llists,[T,3]);

% search list in Cartesian coordinates
list_c = zeros(T,3);
for i = 1 : T
    list_c(i,:) = mtimes(M,list_f(i,:)');
end

% disp(list_f(1,:));
% disp(list_c);

fprintf('\n####################################################################\n')
fprintf('\n Peaks found using Newton-Raphson \n')
% Newton-Raphson scheme
Y_NR = newtonRaphson(list_c, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len, 15, 1e-5);
[~,idx]=unique([round(Y_NR(:,1),4),round(Y_NR(:,2),4),round(Y_NR(:,3),4)],'rows');
peak_NR = Y_NR(idx,:);

TT = length(peak_NR(:,1));
for kk = 1:TT
    fprintf(' Peak %i in FRA = %f %f %f, rho = %f, norm(grad(rho)) = %e  \n',kk, ...
        [peak_NR(kk,1),peak_NR(kk,2),peak_NR(kk,3)], peak_NR(kk,7), peak_NR(kk,8));
end

% prompt = 'Enter the proj2_bench_frac.txt text file: ';
% 
% filename2 = input(prompt,'s');

fprintf('\n####################################################################\n')

fprintf('\n search using Eigenvector following \n')
fprintf('\n length in X direction = %f angstroms\n',len_a)
fprintf('\n length in Y direction = %f angstroms\n',len_b)
fprintf('\n length in Z direction = %f angstroms\n',len_c)

% Grid spacing 0.4 angstrom

% Grid step size in fractional coordinates

GS_a_ef = 0.8*(0.9 - 0.1)/len_a;
GS_b_ef = 0.8*(0.9 - 0.1)/len_b;
GS_c_ef = 0.8*(0.9 - 0.1)/len_c;

% Grid step size (FRA) along CRYSTAL coordinate axes
GS_abc_ef_f = [GS_a_ef, GS_b_ef, GS_c_ef]';
fprintf('\n Grid step size (FRA) along CRYSTAL coordinate axes = \n ')
disp(GS_abc_ef_f');

% Grid step size in cartesian coordinates

GS_abc_ef_c = mtimes(M,GS_abc_ef_f); 
fprintf('\n Grid step size (CAR) along CARTESIAN coordinate axes = \n ')
disp(GS_abc_ef_c');

% All the 3D grid points in Fractional coordinates

F_ef_x = 0.1:GS_a_ef:0.9+GS_a_ef;
F_ef_y = 0.1:GS_b_ef:0.9+GS_b_ef;
F_ef_z = 0.1:GS_c_ef:0.9+GS_c_ef;

len_Fx_ef = length(F_ef_x(:));
len_Fy_ef = length(F_ef_y(:));
len_Fz_ef = length(F_ef_z(:));

fprintf('\n Number of grid points along CRYSTAL coordinate axes = \n')
grid_ef = [len_Fx_ef, len_Fy_ef, len_Fz_ef];
disp(grid_ef);

T_ef = len_Fx_ef*len_Fy_ef*len_Fz_ef; % total length of grid points
fprintf('\n Total number of grid points\n')
disp(T_ef);

F_ef_x = F_ef_x(:);
F_ef_y = F_ef_y(:);
F_ef_z = F_ef_z(:);

% search list in Fractional coordinates
llists_ef = [];
for i = 1:len_Fx_ef
    sublist_ef =[];
    for j = 1:len_Fy_ef
        ssublist_ef=[];
        for k = 1:len_Fz_ef
            ssublist_ef = [ssublist_ef; [F_ef_x(i),F_ef_y(j),F_ef_z(k)]];
        end
        sublist_ef = [sublist_ef; ssublist_ef];
    end
    llists_ef = [llists_ef; sublist_ef];
end
list_ef_f = reshape(llists_ef,[T_ef,3]);

% search list in Cartesian coordinates
list_ef_c = zeros(T_ef,3);
for i = 1 : T_ef
    list_ef_c(i,:) = mtimes(M,list_ef_f(i,:)');
end


fprintf('\n####################################################################\n')
fprintf('\n Peaks found using Eigenvector following \n')
% Eigenvector following scheme
Y_EF = EigVectorFollow(list_ef_c, A_hkl, B_hkl, H, K, L, Vc, M_inv, data_len, 15, 1e-5);
[~,idx]=unique([round(Y_EF(:,1),4),round(Y_EF(:,2),4),round(Y_EF(:,3),4)],'rows');
peak_EF = Y_EF(idx,:);

TT = length(peak_EF(:,1));
for kk = 1:TT
    fprintf(' Peak %i in FRA = %f %f %f, rho = %f, norm(grad(rho)) = %e  \n',kk, ...
        [peak_EF(kk,1),peak_EF(kk,2),peak_EF(kk,3)], peak_EF(kk,7), peak_EF(kk,8));
end

fprintf('\n####################################################################\n')
toc
%%%%%%%%%%%%%%%%%%% functions used in this project %%%%%%%%%%%%%%%%%%%%%%%%

function rhof = rhoxf(x, A, B, h, k, l, V, len)
Sumf = 0;
for i = 1:len
    Sumf = Sumf + A(i)*cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)))...
    + B(i)*sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
    rhof = (1/V)*2*Sumf;
end
end

function rhoc = rhoxc(xx, A, B, h, k, l, V, Minv, len)
Sumc = 0;
x = mtimes(Minv,xx);
for i = 1:len
    Sumc = Sumc + A(i)*cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)))... 
    + B(i)*sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
    rhoc = (1/V)*2*Sumc;
end
end

function gradd = gradient_c(xx, A, B, h, k, l, V, Minv, len)
gradd = zeros(3,1);
SumG = zeros(3,1);
x = mtimes(Minv,xx);
for j =1:length(SumG)
    for i = 1:len
        SumG(j) = SumG(j) + ...
            2*pi*(Minv(1,j).*h(i) + Minv(2,j).*k(i) + Minv(3,j).*l(i))*B(i)...
            *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
            2*pi*(Minv(1,j).*h(i) + Minv(2,j).*k(i) + Minv(3,j).*l(i))*A(i)...
            *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));      
    end
    gradd(j) = (1/V)*2*SumG(j);
end
end


function hess_c = hessian_c(xx, A, B, h, k, l, V, Minv, len)
hess_c = zeros(3);
x = mtimes(Minv,xx);

sum11 = 0;
for i = 1:len
    sum11 = sum11 + ...
        -4*pi^2*A(i)*(Minv(1,1).*h(i) + Minv(2,1).*k(i) + Minv(3,1).*l(i))^2 ...
        *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
        4*pi^2*B(i)*(Minv(1,1).*h(i) + Minv(2,1).*k(i) + Minv(3,1).*l(i))^2 ...
        *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
end

sum22 = 0;
for i = 1:len
    sum22 = sum22 + ...
        -4*pi^2*A(i)*(Minv(1,2).*h(i) + Minv(2,2).*k(i) + Minv(3,2).*l(i))^2 ...
        *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
        4*pi^2*B(i)*(Minv(1,2).*h(i) + Minv(2,2).*k(i) + Minv(3,2).*l(i))^2 ...
        *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
end

sum33 = 0;
for i = 1:len
    sum33 = sum33 + ...
        -4*pi^2*A(i)*(Minv(1,3).*h(i) + Minv(2,3).*k(i) + Minv(3,3).*l(i))^2 ...
        *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
        4*pi^2*B(i)*(Minv(1,3).*h(i) + Minv(2,3).*k(i) + Minv(3,3).*l(i))^2 ...
        *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
end

sum12 = 0;
for i = 1:len
    sum12 = sum12 + ...
        -4*pi^2*A(i)*(Minv(1,1).*h(i) + Minv(2,1).*k(i) + Minv(3,1).*l(i)) ...
        *(Minv(1,2).*h(i) + Minv(2,2).*k(i) + Minv(3,2).*l(i)) ...
        *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
        4*pi^2*B(i)*(Minv(1,1).*h(i) + Minv(2,1).*k(i) + Minv(3,1).*l(i)) ...
        *(Minv(1,2).*h(i) + Minv(2,2).*k(i) + Minv(3,2).*l(i)) ...
        *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
end

sum13 = 0;
for i = 1:len
    sum13 = sum13 + ...
        -4*pi^2*A(i)*(Minv(1,1).*h(i) + Minv(2,1).*k(i) + Minv(3,1).*l(i)) ...
        *(Minv(1,3).*h(i) + Minv(2,3).*k(i) + Minv(3,3).*l(i)) ...
        *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
        4*pi^2*B(i)*(Minv(1,1).*h(i) + Minv(2,1).*k(i) + Minv(3,1).*l(i)) ...
        *(Minv(1,3).*h(i) + Minv(2,3).*k(i) + Minv(3,3).*l(i)) ...
        *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
end

sum23 = 0;
for i = 1:len
    sum23 = sum23 + ...
        -4*pi^2*A(i)*(Minv(1,2).*h(i) + Minv(2,2).*k(i) + Minv(3,2).*l(i)) ...
        *(Minv(1,3).*h(i) + Minv(2,3).*k(i) + Minv(3,3).*l(i)) ...
        *cos(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i))) - ...
        4*pi^2*B(i)*(Minv(1,2).*h(i) + Minv(2,2).*k(i) + Minv(3,2).*l(i)) ...
        *(Minv(1,3).*h(i) + Minv(2,3).*k(i) + Minv(3,3).*l(i)) ...
        *sin(2*pi*(x(1)*h(i) + x(2)*k(i) + x(3)*l(i)));
end
hess_c(1,1) = (1/V)*2*sum11;
hess_c(1,2) = (1/V)*2*sum12;
hess_c(1,3) = (1/V)*2*sum13;
hess_c(2,1) = (1/V)*2*sum12;
hess_c(2,2) = (1/V)*2*sum22;
hess_c(2,3) = (1/V)*2*sum23;
hess_c(3,1) = (1/V)*2*sum13;
hess_c(3,2) = (1/V)*2*sum23;
hess_c(3,3) = (1/V)*2*sum33;
end

function [eigVe, eigVa] = eigg(xx, A, B, h, k, l, V, Minv, len)
hess_c = hessian_c(xx, A, B, h, k, l, V, Minv, len);
[eigVe, eigVa] = eig(hess_c);
end

function Y = newtonRaphson(xx, A, B, h, k, l, V, Minv, len, N, tol)
Y = [];

peak_ff = [];
eigVa_peak_ff = [];
rho_peak_ff = [];
n_grad_peak_ff = [];

last_peak_ff = [];
last_eigVa_peak_ff = [];
last_rho_peak_ff = [];
last_n_grad_peak_ff = [];

x = xx';

for i = 1:length(xx)
    xcc = [];
    xff = [];
    last_xff = [];
    last_rho_c = [];
    last_n_grad_c = [];
    last_eigVa_c = [];
    rho_cc = [];
    n_grad_cc = [];
    eigVa_cc = [];
    xc_old = x(:,i);
    xf_old = mtimes(Minv,xc_old);
    
    while length(xcc)<= N
        xcc = [xcc; xc_old];
        xff = [xff; xf_old];
        grad_c = gradient_c(xc_old, A, B, h, k, l, V, Minv, len);
        hess_c = hessian_c(xc_old, A, B, h, k, l, V, Minv, len);
        rho_c = rhoxc(xc_old, A, B, h, k, l, V, Minv, len);
        rhoff = rhoxf(xf_old, A, B, h, k, l, V, len); % new edit       
        n_grad_c = norm(grad_c);
        [eigVe, eigVa] = eigg(xc_old, A, B, h, k, l, V, Minv, len);
        eigVa_c = diag(eigVa);
        
        eigVa_cc = [eigVa_cc; eigVa_c];
        rho_cc = [rho_cc; rho_c];
        n_grad_cc = [n_grad_cc; n_grad_c];
        
        hh = linsolve(-hess_c ,grad_c);
        if n_grad_c < tol && rho_c < 2
            break
        end 
        
        if norm(hh) > 0.25
            n_hh = norm(hh);
            hh = (0.25/n_hh)*hh;
        end
        
        xc_new = xc_old + hh;
        xf_new = mtimes(Minv,xc_new);
        
        xc_old = xc_new;
        xf_old = xf_new;        
    end            
r_xcc = reshape(xcc,[],1);
r_xff = reshape(xff,[],1);

r_eigVa_cc = reshape(eigVa_cc,[],1);
r_rho_cc = reshape(rho_cc,[],1);
n_grad_cc = reshape(n_grad_cc,[],1);

peak_ff = [peak_ff;r_xff];
eigVa_peak_ff = [eigVa_peak_ff;r_eigVa_cc];
rho_peak_ff = [rho_peak_ff;r_rho_cc];
n_grad_peak_ff = [n_grad_peak_ff;n_grad_cc];

last_xff_a = reshape(peak_ff,3,[]);
last_xff = [last_xff_a(1,end), last_xff_a(2,end), last_xff_a(3,end)];
last_eigVa_c_a = reshape(eigVa_peak_ff,3,[]);
last_eigVa_c = [last_eigVa_c_a(1,end), last_eigVa_c_a(2,end), last_eigVa_c_a(3,end)];
last_rho_c_a = reshape(rho_peak_ff,1,[]);
last_rho_c = last_rho_c_a(1,end);
last_n_grad_c_a = reshape(n_grad_peak_ff,1,[]);
last_n_grad_c = last_n_grad_c_a(1,end);

last_peak_ff = [last_peak_ff;last_xff];
last_eigVa_peak_ff = [last_eigVa_peak_ff;last_eigVa_c];
last_rho_peak_ff = [last_rho_peak_ff;last_rho_c];
last_n_grad_peak_ff = [last_n_grad_peak_ff;last_n_grad_c];
end

YY = [last_peak_ff last_eigVa_peak_ff last_rho_peak_ff last_n_grad_peak_ff]; 

for k = 1:length(last_peak_ff)
    YYY = [YY(k,1),YY(k,2),YY(k,3),YY(k,4),YY(k,5),YY(k,6),YY(k,7),YY(k,8)];
    if all([YY(k,4), YY(k,5), YY(k,6)] < 0) && YY(k,7) > 2
        Y = [Y;YYY];
    end
end

end

function Y_ef = EigVectorFollow(xx, A, B, h, k, l, V, Minv, len, N, tol)
Y_ef = [];

peak_ff = [];
eigVa_peak_ff = [];
rho_peak_ff = [];
n_grad_peak_ff = [];

last_peak_ff = [];
last_eigVa_peak_ff = [];
last_rho_peak_ff = [];
last_n_grad_peak_ff = [];

x = xx';
for i = 1:length(xx)
    B_prime = zeros(4);
    h_prime = zeros(4,1);
    xcc = [];
    xff = [];
    last_xff = [];
    last_rho_c = [];
    last_n_grad_c = [];
    last_eigVa_c = [];
    rho_cc = [];
    n_grad_cc = [];
    eigVa_cc = [];
    xc_old = x(:,i);
    xf_old = mtimes(Minv,xc_old);
    
    while length(xcc)<= N
        xcc = [xcc; xc_old];
        xff = [xff; xf_old];
        grad_c = gradient_c(xc_old, A, B, h, k, l, V, Minv, len);
        hess_c = hessian_c(xc_old, A, B, h, k, l, V, Minv, len);
        rho_c = rhoxc(xc_old, A, B, h, k, l, V, Minv, len);
        n_grad_c = norm(grad_c);
        [eigVe, eigVa] = eigg(xc_old, A, B, h, k, l, V, Minv, len);
        eigVa_c = diag(eigVa);
        eigVe_1c = eigVe(:,1);
        eigVe_2c = eigVe(:,2);
        eigVe_3c = eigVe(:,3);
        F1 = mtimes(eigVe_1c',grad_c);
        F2 = mtimes(eigVe_2c',grad_c);
        F3 = mtimes(eigVe_3c',grad_c);
        
        bb  = sort(eigVa_c);
        
        B_prime(1,1) = bb(1);  
        B_prime(2,2) = bb(2);
        B_prime(3,3) = bb(3);
        B_prime(1,4) = F1;
        B_prime(2,4) = F2;
        B_prime(3,4) = F3;
        B_prime(4,1) = F1;
        B_prime(4,2) = F2;
        B_prime(4,3) = F3;
        
        hh = linsolve(-hess_c ,grad_c);
        
        h_prime(1,1) = hh(1);
        h_prime(2,1) = hh(2);
        h_prime(3,1) = hh(3);
        h_prime(4,1) = 1;
        
        [lamVe, lamVa] = eig(B_prime);
        
        eig_max = max(diag(lamVa));
        
        eigVa_cc = [eigVa_cc; eigVa_c];
        rho_cc = [rho_cc; rho_c];
        n_grad_cc = [n_grad_cc; n_grad_c];
        
        hh = -((mtimes(F1,eigVe_1c)/(eigVa_c(1)-eig_max)) + ...
            (mtimes(F2,eigVe_2c)/(eigVa_c(2)-eig_max)) + ...
            (mtimes(F3,eigVe_3c)/(eigVa_c(3)-eig_max)));
        
        
        if n_grad_c < tol && rho_c < 2.5
            break
        end 
        
        if norm(hh) > 0.25
            n_hh = norm(hh);
            hh = (0.25/n_hh)*hh;
        end
        
        xc_new = xc_old + hh;
        xf_new = mtimes(Minv,xc_new);
        
        xc_old = xc_new;
        xf_old = xf_new;        
    end            
r_xcc = reshape(xcc,[],1);
r_xff = reshape(xff,[],1);

r_eigVa_cc = reshape(eigVa_cc,[],1);
r_rho_cc = reshape(rho_cc,[],1);
n_grad_cc = reshape(n_grad_cc,[],1);

peak_ff = [peak_ff;r_xff];
eigVa_peak_ff = [eigVa_peak_ff;r_eigVa_cc];
rho_peak_ff = [rho_peak_ff;r_rho_cc];
n_grad_peak_ff = [n_grad_peak_ff;n_grad_cc];

last_xff_a = reshape(peak_ff,3,[]);
last_xff = [last_xff_a(1,end), last_xff_a(2,end), last_xff_a(3,end)];
last_eigVa_c_a = reshape(eigVa_peak_ff,3,[]);
last_eigVa_c = [last_eigVa_c_a(1,end), last_eigVa_c_a(2,end), last_eigVa_c_a(3,end)];
last_rho_c_a = reshape(rho_peak_ff,1,[]);
last_rho_c = last_rho_c_a(1,end);
last_n_grad_c_a = reshape(n_grad_peak_ff,1,[]);
last_n_grad_c = last_n_grad_c_a(1,end);

last_peak_ff = [last_peak_ff;last_xff];
last_eigVa_peak_ff = [last_eigVa_peak_ff;last_eigVa_c];
last_rho_peak_ff = [last_rho_peak_ff;last_rho_c];
last_n_grad_peak_ff = [last_n_grad_peak_ff;last_n_grad_c];
end

YY = [last_peak_ff last_eigVa_peak_ff last_rho_peak_ff last_n_grad_peak_ff]; 

for k = 1:length(last_peak_ff)
    YYY = [YY(k,1),YY(k,2),YY(k,3),YY(k,4),YY(k,5),YY(k,6),YY(k,7),YY(k,8)];
    if all([YY(k,4), YY(k,5), YY(k,6)] < 0) && YY(k,7) > 2.5 && YY(k,8) < 1e-5 
        Y_ef = [Y_ef;YYY];
    end
end

end

