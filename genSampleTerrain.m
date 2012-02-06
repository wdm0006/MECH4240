function out=genSampleTerrain(numrows, numcols)

out=zeros(numrows, numcols);
north_gradient=.0015*(pi/180);
east_gradient=.001*(pi/180);
for j=2:numrows
   out(j,:)=out(j-1,:)+sin(north_gradient); %-(rand(1,1)/10);
end
    
for j=2:numcols
   out(:,j)=out(:,j-1)+sin(east_gradient); %-(rand(1,1)/10);
end
 

    