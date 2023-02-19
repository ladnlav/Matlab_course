clear all; clc; close all;

%% Analisys of text from file
% to-do
% 1) Reading the file as a text of array of chars
% 2) Create array of cells which consist of three columns
% "char"->"amount of meetings"->"probobilities"
% 3) Save the chars and probobalities to file *.mat and *.xls as the cell
% variables. Name the files shouold be:
% Data_Analisys.mat
% Data_Analisys.xls
% Data_Analisys.png
% 4) Plot the distribution of probability of symbols in text. 
% Be careful to the labels on the axis.
% Recommendation use xticks(), xticklabels().
% 5) Save the plot as figure and PNG image with resolution at least 400 px. The name
% of files should be: Data_Analisys.png

%% Reading the file
% TO-DO 1
% Read the file from *.txt as a char stream
fid = fopen('Textvar6.txt', 'r+');
[Char_from_File, Size_from_file] = fscanf(fid, '%c'); %чтение текстового ф.
fclose(fid);

%% Analysis
% TO-DO 2 
% Use ony char from file
% Use  lowercase string
% Try to use the "Cell" as a data containers;
% Name the varible Data_Analisys
% The cell should consist of 3 columns:
% "Symbol" | "Amount of meeting" | "Probolitie"

% You can use only 1 cycle for this task
% Avoid the memmory allocation in cycle

str=lower(Char_from_File);
ustr=unique(str);
Data_Analisys = cell(strlength(ustr),3);
itter=1;

for ch = ustr
    Data_Analisys{itter,1}=ch;
    Data_Analisys{itter,2}=count(str,ch);
    Data_Analisys{itter,3}=Data_Analisys{itter,2}/strlength(str);

    itter=itter+1;
end

%% Plot Data
% TO-DO 3
% Illustrate the results from Analysis block
% There should be lable of axis, title, grid

bar([Data_Analisys{:,3}])
X_l=Data_Analisys(:,1);
X_l = regexprep(X_l, '\n', sprintf('newline'));
X_l = regexprep(X_l, '\r', sprintf('begline'));
X_l = regexprep(X_l, '\ ', sprintf('space'));
xticks(1:size(Data_Analisys,1))
xticklabels(X_l)
grid on;
title("The distribution of probability of symbols in text");
xlabel('Symbols');
ylabel('Probability');


%% Save the file
% TO-DO 4
% Save the figure as Data_Analisys.fig
saveas(gcf,'Data_Analisys.fig')
% Save the figure as image Data_Analisys.png
saveas(gcf,'Data_Analisys.png')
% Save the data as MAT-file Data_Analisys.mat
save('Data_Analisys.mat',"Data_Analisys")
% Save the data as Excel table Data_Analisys.xls
writecell(Data_Analisys, 'Data_Analisys.xls');
%% Discrete lvls of probability

sorted_DA=sort([Data_Analisys{:,3}]);
bar(sorted_DA)
xticks(1:size(Data_Analisys,1))
grid on;
title("The Discret lvls of probability");
xlabel('Symbols');
ylabel('Probability');


