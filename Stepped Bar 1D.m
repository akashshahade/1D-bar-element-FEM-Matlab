
%-------------------------------------------------------------------%
% Software by -  Akash S. Shahade
% Title - Program for FE Analysis of 1D Bar.
%-------------------------------------------------------------------%
clear;
clc;
fprintf(' \nFEM ANALYSIS OF 1-D BAR WITH DIFFERENT CROSS-SECTIONS. \n\n');

%------------------------------------------------------%
% STEP 01 - PRE-PROCESSING
%------------------------------------------------------%


ele_nod=[1 2;2 3];                % Elements are connected with these nodes
nod_coor=[0 0;0 0.5;0 1];         % Node coordinates
num_ele=size(ele_nod,1);          % Number of elements
ele_dof=[1 2;2 3];                % D.O.F  associated with Nodes [1 2;2 3]
num_nod=3;                        % Number of Nodes
dof = 1;                          % D.O.F per node

displacement = zeros(dof*num_nod,1);    % Zero Matrix for Displacement
force = zeros(dof*num_nod,1);           % Zero Matrix for Force
stiffness = zeros(dof*num_nod);         % Zero Matrix for Stiffness

A(1)=200;                   % Area of Element 01 (mm^2)
A(2)=180;                   % Area of Element 02 (mm^2)
L(1)=20;                    % Length of Element 01 (mm)
L(2)=60;                    % Length of Element 02 (mm)
E(1)=200*10^3;              % Young's Modulus of Element 01 (Mpa)
E(2)=120*10^3;              % Young's Modulus of Element 02 (Mpa)

%------------------------------------------------------%
% Stiffness matrix calculation & ASSEMBLY
%------------------------------------------------------%


for e=1:num_ele                               % For 1 to Number of elements
    
    k = ((A(e)*E(e))/L(e))*[1 -1;-1 1];

   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);

    for i=1:2
        for j=1:2
                                              % Assembly of Global Matrix
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);

        end
    end
end


force(3)=10000;

fprintf('Global Stiffness Matrix: \n');
disp(stiffness);
fprintf('\n Global Load Vector: \n');
disp(force);
fprintf('\n------------------------------------------------\n');


%------------------------------------------------------%
% Boundary Conditions
%------------------------------------------------------%

fixed_dof = 1;                       % Constrained D.O.F.
k=stiffness;
k(fixed_dof,:)=[];                   % Eliminating Rows
k(:,fixed_dof)=[];                   % Eliminating Columns

f=force;
f(fixed_dof,:)=[];                   % Eliminating Rows

%------------------------------------------------------%
% STEP 02 - SOLVE
%------------------------------------------------------%

q = k\f ;

%------------------------------------------------------%
% STEP 03 - POST-PROCESSING
%------------------------------------------------------%

displacement=[0;q];                         % Displacement Vector

reaction = stiffness*displacement - force;  % Calculate Reaction Forces

Node = [1;2;3];
qxy = {'u1';'u2';'u3'};
Q_mm = displacement;
T=table(Node,qxy,Q_mm);
disp(T);

fprintf('\n------------------------------------------------\n');

F = {'R1x';'R2x';'R3x'};
Reaction_N = reaction;
J = table(F,Reaction_N);
disp(J);
fprintf('END OF PROGRAM.\n');
