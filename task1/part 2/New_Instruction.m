function L = New_Instruction(Mat)

%first cycle
Mat(26:end, :) = abs(Mat(26:end, :));
%second cycle
col_idx = mod(1:size(Mat,2), 9) ~= 0;
Mat(:, col_idx) = 7 * Mat(:, col_idx);
%third cycle
%also possible:
%[i,j] = find(Mat<-50);
%Mat(i,j)=-Mat(i,j);
Mat(Mat<-50)=-Mat(Mat<-50);

L=Mat;
end