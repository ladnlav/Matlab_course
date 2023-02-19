% Clear workspace. Clear command window. Close all figure
clear all; clc; close all;
%% task
% 0) Create a function
% 1) Create the random matrix
% 2) Analyse the code. Insert the calculation of runtime of code
% 3) rewrite the code in more optimised form in matlab
% 4) Provide the evidence that results matrix and legacy matrix is the same
% 5) calculate the runtime of new code. Compare it with legacy code. Make
% an conclusion about code. Which one is the more optimised? Which code do
% you suggest to use in matlab? And why?
%% Config the system

% Fixed random generator
rng(110);
% TO-DO 1%
% Create function, which generate 
% Input_Matrix matrix 200-to-600 size and
% with normal distributed numbers from -100 to 22
%    

left=-100;
right=22;
rows=200;
col=600;

Input_Matrix = Matrix_generator(left,right,rows,col);  
Legacy_Matrix = Input_Matrix;
Ethalon_Matrix = Input_Matrix;
%% Run legacy code
% TO-DO 2
% Measure the runtime of current function
tic
Legacy_output_Matrix = Legacy_Instruction(Input_Matrix);

% Save the runtime in variable
Time_legacy_code = toc;

%% Run optimised code
% TO-DO 3
% Measure the runtime of your function
% Create function New_Instruction()
% Rewrite and optimised function Legacy_Instruction()
% Use matrix operation, logical indexing
% Try not to use the cycle
tic
Optimised_Output_Matrix = New_Instruction(Input_Matrix);
% Save the runtime in variable
Time_Optimised_code = toc;
    
%% Checking the work of student
% TO-DO 4
% Compare the matrix and elapsed time for instruction
% Matrix must be equal each other, but the runtime sill be different

%Optimised_Output_Matrix(1,1)=0;

if isequal(Optimised_Output_Matrix,Legacy_output_Matrix)
    disp('Matrices equal!')
else
    disp('Smth WRONG!')
end

% Runtime comparison
if Time_legacy_code-Time_Optimised_code==0
    disp('Smth WRONG!')
else
    disp(['Optimized code is faster in ' num2str((Time_legacy_code/Time_Optimised_code)) ' times!']);
end

% Comparison of matrix
    % Matrix size and value 
%% Statistical analysis
% Determine the number of function launches - n
n=100000;
Input_Matrix=Matrix_generator(left,right,rows,col);

% Create arrays to store the execution time of functions
t1 = zeros(n, 1);
t2 = zeros(n, 1);

% Run the functions several times
for i = 1:n

    tic;
    Legacy_output_Matrix = Legacy_Instruction(Input_Matrix);
    t1(i) = toc;
    
    tic;
    Optimised_Output_Matrix = New_Instruction(Input_Matrix);
    t2(i) = toc;
end

ratio=t1./t2;
% Output the results

plot(1:n, ratio, 'bo-'); 
xlabel('Number of the launch');
ylabel('Ratio of times');
title('Dependence of the ratio of function execution times on the number of launches');
disp(['(Average) Optimized code is faster in ' num2str(mean(ratio)) ' times!']);
disp(['Standard deviation of the execution time of the ratio: ', num2str(std(ratio))]);


%% Function discribing

function Output_Matrix = Legacy_Instruction(Matrix)
   
    for itter_rows = 1 : size(Matrix,1)
        for itter_column = 1 : size(Matrix,2)
            if itter_rows > 25
                Matrix(itter_rows,itter_column) = abs(Matrix(itter_rows,itter_column));
            end
        end
    end

   for itter_rows = 1 : size(Matrix,1)
        for itter_column = 1 : size(Matrix,2)
            if mod(itter_column,9) ~= 0
                Matrix(itter_rows,itter_column) = 7*Matrix(itter_rows,itter_column);
            end
        end
   end

   for itter_rows = 1 : size(Matrix,1)
        for itter_column = 1 : size(Matrix,2)
            if Matrix(itter_rows,itter_column) < -50
                Matrix(itter_rows,itter_column) = -Matrix(itter_rows,itter_column);
            end
        end
    end

    Output_Matrix = Matrix;
end