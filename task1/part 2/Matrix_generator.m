function L = Matrix_generator(left,right,rows,col)
L = (right-left).*rand(rows, col)+ left;
end