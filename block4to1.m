%% block8to4: function description
function result = block4to1(image1,image2,image3,image4)

[W,L] = size(image1);
s = 8*ones(1,W/8);

imagecell1 = mat2cell(image1,s,s);
imagecell2 = mat2cell(image2,s,s);
imagecell3 = mat2cell(image3,s,s);
imagecell4 = mat2cell(image4,s,s);

tempcell = cell(size(imagecell1)*2);

[W,L] = size(imagecell1);

for i = 1:W
	for j = 1:L
		tempcell{2*(i-1)+1,2*(j-1)+1} = imagecell1{i,j};
		tempcell{2*(i-1)+1,2*(j-1)+2} = imagecell2{i,j};
		tempcell{2*(i-1)+2,2*(j-1)+2} = imagecell3{i,j};
		tempcell{2*(i-1)+2,2*(j-1)+1} = imagecell4{i,j};
	end
end
result = cell2mat(tempcell);

end

