function pic_normalized = normalize(pic)

%normalize picture: scale from -1 to 1
min_value = min(min(min(pic)));
max_value = max(max(max(pic)));
middle = (max_value+min_value)/2;
pic_normalized = (pic-middle)/(max_value-middle);

end