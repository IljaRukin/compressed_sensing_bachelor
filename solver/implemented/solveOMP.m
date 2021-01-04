function x = OMP(pic_sparse,dim,indices,sparsity,start_lsqlin_OptimalityTolerance,final_lsqlin_OptimalityTolerance,accuracy,AA,AT)

l2sqr = @(A) sum(sum(sum(A.*conj(A))));
l2_samples = l2sqr(pic_sparse);

result_positions = [];
result_values = [];

residue = pic_sparse;
A=[];
null_vec=zeros(dim);

%options for !!!coarse!!! lsqlin coefficient value calculation
tolerance_steps = 11;
tolerances = logspace(log10(start_lsqlin_OptimalityTolerance),log10(final_lsqlin_OptimalityTolerance),tolerance_steps);
change_tolerance = 1;
opts = optimoptions('lsqlin', 'OptimalityTolerance', start_lsqlin_OptimalityTolerance);

nr_coeff=0;
error=1;
while nr_coeff<sparsity && error>accuracy
    nr_coeff = nr_coeff+1;
    
    %calculate correlation
    expanded_vec = zeros(dim);
    expanded_vec(indices) = residue;
    correlation = AA(expanded_vec);
    
    %find entry with highest correlation
    [~,max_position] = max(reshape(abs(correlation),prod(dim),1));
    max_value = correlation(max_position) * prod(dim)/numel(indices);
    
    %test if that member is already taken into account
    if ismember(max_position,result_positions)
        %improve accuracy and reexecute loop
        change_tolerance = change_tolerance +1;
        if change_tolerance > tolerance_steps
            fprintf('-----stop: no change in accuracy by adding coefficients----- \n');
            fprintf('-----lsqlin optimality tolerance limit reached !----- \n');
            break
        end
        opts = optimoptions('lsqlin', 'OptimalityTolerance', tolerances(change_tolerance));
        continue
    end
    
    %add new coefficient
    result_positions = [result_positions,max_position];
    
    %setup equation for accounted coefficients and solve with lsqlin
    null_vec(max_position)=1;
    temp = AT(null_vec);
    A=[A,temp(indices)];
    null_vec=zeros(dim);
    %A=reverse_basis_sparse(:,result_positions);
    result_values = lsqlin(A,pic_sparse,[],[],[],[],[],[],[result_values;max_value],opts);
    residue = pic_sparse - A*result_values;
    error = l2sqr(residue)/l2_samples;
end

%final lsqlin coefficient value calculation (better accuracy)
opts = optimoptions('lsqlin', 'OptimalityTolerance', final_lsqlin_OptimalityTolerance);
result_values = lsqlin(A,pic_sparse,[],[],[],[],[],[],result_values,opts);
residue = pic_sparse - A*result_values;

%calculate error
error = l2sqr(residue)/l2_samples;

%reconstruct coefficient matrix
x=zeros(dim);
x(result_positions)=result_values;

fprintf('-----termination----- \n');
fprintf('sparsity = %i \n',numel(result_positions));
fprintf('error = %i \n',error);
fprintf('-------------------- \n');

end