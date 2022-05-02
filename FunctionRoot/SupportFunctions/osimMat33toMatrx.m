function Mat = osimMat33toMatrx(Mat33)

    for i = 1:Mat33.nrow()
        for j = 1:Mat33.ncol()
            Mat(i, j) = Mat33.get(i-1, j-1);
        end
    end
    
end