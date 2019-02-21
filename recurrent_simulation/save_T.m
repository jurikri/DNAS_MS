T = zeros(60,60);
T(:,:) = 255;
T(15:20,20:40) = 12;
T(20:end-10,26:33) = 12;
T = uint8(T);

for row = 1:size(T,1)
    T2(row,:) = T(size(T,1)+1-row,:);
end

imshow(T2)
imwrite(T2,'t.bmp')