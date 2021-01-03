function pic_normalized = normalization(pic)

%normalize picture: mean=0 & variance=1
mean = sum(sum(sum(pic)))/numel(pic);
variance = sum(sum(sum((pic-mean).^2)))/numel(pic);
pic_normalized = (pic - mean)/sqrt(variance);

end