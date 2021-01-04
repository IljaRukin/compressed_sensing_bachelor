function x = IHT(pic_sparse,dim,indices,accuracy,min_change,maxiter,min_step,step_search,start_sparsity,max_sparsity,sparsity_step,AA,AT)

%IHT approximations for sparsity, step and stopping criteria
%
%###stop
%if error smaller than accuracy
%step smaller than min_step
%
%###increase sparsity
%if error does not change by min_change
%if step is smaller than min_step an sparsity greater or equal than max_sparsity
%
%###new step
%error after gradient step smaller than before => use this step (latch1)
%continue decreasing step (update best_opt) untill error increases again or does not change greater than min_change (latch2)


l2sqr = @(A) sum(sum(sum(A.*conj(A))));
l2_samples = l2sqr(pic_sparse);

%setup termination measures
iter = 1;
step = 1;%*10
coeff = start_sparsity;

%first guess
expanded = zeros(dim);
expanded(indices) = pic_sparse;
x_next = AT(expanded);
[~,maxindices] = sort(reshape(abs(x_next),prod(dim),1),'descend');
x = zeros(dim);
x(maxindices(1:coeff)) = x_next(maxindices(1:coeff));

%for line search set starting optima
temp = AT(x);
error = l2sqr(temp(indices)-pic_sparse)/l2_samples;
opt_temp = error;

while iter<maxiter
    
    %Calculation
    trafo = AT(x);
    difference = pic_sparse - trafo(indices);
    
    %test accuracy -> early termination
    last_error = error;
    error = l2sqr(difference)/l2_samples;
    if error<accuracy
        fprintf('-----stop: accuracy reached----- \n');
        break
    end
    
    %if error does not change by min_change => increase sparsity
    %and if max_sparsity reached => stop
    if last_error < (1+min_change)*error && coeff < max_sparsity && iter~=1
        coeff = coeff + sparsity_step;
    elseif last_error < (1+min_change)*error && coeff >= max_sparsity
        fprintf('-----stop: max sparsity reached and error change smaller than min_change----- \n');
        break
    end
    
    %resume Calculation
    expanded = zeros(dim);
    expanded(indices) = difference;
    gradient = AA(expanded);
    l2_change = norm(gradient);
    
    %approximate line search for convex optima
    step = step*10;
    latch1 = 0;
    latch2 = 0;
    %set optima for current x, it is the optima to beat
    best_opt = error; %error to beat !!!
    while true
        previous_opt = opt_temp;
        if step < min_step && coeff < max_sparsity %accuraccy cant be improve => increase sparsity
            coeff = coeff + sparsity_step;
			step = 1;
        elseif step < min_step %step too small && maximum sparsity reached => terminate
            break
        end
        x_temp = x + step/l2_change*gradient;
        [~,maxindices] = sort(reshape(abs(x_temp),prod(dim),1),'descend');
        temp_x = zeros(dim);
        temp_x(maxindices(1:coeff)) = x_temp(maxindices(1:coeff));
        temp = AT(temp_x);
        opt_temp = l2sqr(temp(indices)-pic_sparse)/l2_samples;
        if opt_temp>previous_opt %optima starts to increase again => prepare for step
            latch1 = 1;
        end
        if opt_temp < best_opt %better optima found => prepare for step
            best_opt = opt_temp;
            best_step = step;
            latch2 = 1;
        elseif opt_temp*(1+min_change) > previous_opt %optima does not change much => prepare for step
            latch1 = 1;
        end
        if latch1 && latch2 %if all latches set => take step
            break
        end
        step = step*step_search;
    end
    
    if step < min_step && coeff < max_sparsity
        fprintf('-----stop: step too small----- \n');
        break
    end
    
    x_next = x + best_step/l2_change*gradient;
    
    %Hard Thresholding
    x = zeros(dim);
    [~,maxindices] = sort(reshape(abs(x_next),prod(dim),1),'descend');
    x(maxindices(1:coeff)) = x_next(maxindices(1:coeff));
    
    iter = iter +1;
    
end
if iter>=maxiter
    fprintf('-----stop: iteration limit reached----- \n');
end
fprintf('-----termination----- \n');
fprintf('iterations = %i \n',iter);
fprintf('sparsity = %i \n',coeff);
fprintf('error = %i \n',error);
fprintf('-------------------- \n');

end