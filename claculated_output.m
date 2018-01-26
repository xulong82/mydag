load sachsDiscretizedData;


a=[1000 500 250 100 50];
output={};

parfor i=1:5
data1=data(:,1:a(i));
output{i}=Score_all(data1,ns);
end


