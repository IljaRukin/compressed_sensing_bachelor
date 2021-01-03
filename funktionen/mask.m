function mask = mask(dim,percent_ones)

%create sampling mask

pixels = prod(dim);
samples = round(pixels*percent_ones);

seed=1;
rng('default');
rng(seed);

mask = zeros(dim);
mask(randperm(numel(mask),samples))=1;
