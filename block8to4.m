%% block8to4: function description
function [result1,result2,result3,result4] = block8to4(image)

[W,L] = size(image);
s = 8*ones(1,W/8);

imagecell = mat2cell(image,s,s);
[tempcell1,tempcell2,tempcell3,tempcell4] = deal(cell(size(imagecell)/2));

[W,L] = size(tempcell4);
for i = 1:W
	for j = 1:L
		tempcell1{i,j} = imagecell{2*(i-1)+1,2*(j-1)+1};
		tempcell2{i,j} = imagecell{2*(i-1)+1,2*(j-1)+2};
		tempcell3{i,j} = imagecell{2*(i-1)+2,2*(j-1)+2};
		tempcell4{i,j} = imagecell{2*(i-1)+2,2*(j-1)+1};
	end
end
result1 = cell2mat(tempcell1);
result2 = cell2mat(tempcell2);
result3 = cell2mat(tempcell3);
result4 = cell2mat(tempcell4);
end

