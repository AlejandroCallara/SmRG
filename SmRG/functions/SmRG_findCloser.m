function [sample_values,found,stored_init_cond,proposed]=...
    find_closer(actual_values,sample_values,stored_init_cond,...
    proposed,threshold)
%input
%actual_values  a vector containing the

if ~isrow(actual_values)
    actual_values = actual_values';
end

if isempty(stored_init_cond)
    sample_values=[];%actual_values;
    stored_init_cond=proposed;
    found=0;
else
    Num_Samples=size(sample_values,1)
    distances=zeros(Num_Samples,1);
    for kk=1:Num_Samples
        %find the common x axis
        min_c=min(min(sample_values(kk,:)),min(actual_values));
        max_c=max(max(sample_values(kk,:)),max(actual_values));
        num_bins=floor(log2(length(sample_values))+1);
        bin_limits=linspace(min_c,max_c,num_bins);
        p1=histc(sample_values(kk,:),bin_limits);
        p2=histc(actual_values,bin_limits);
        d1=abs(KLDiv(p1,p2));
        distances(kk)=d1;
    end
    [value_min,ind_min]=min(distances);
    if value_min<threshold
        initial_values=stored_init_cond(ind_min,:);
        found=1;
        proposed = initial_values;
    else
        found=0;
        proposed = [];
    end
    
end






