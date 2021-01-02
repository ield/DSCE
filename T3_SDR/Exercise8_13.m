%% Exercise 8_13
% Find the position of the marker in the data

% Two posible marks descripted in the chapter
% mark = ones(1, 7);
mark = [1 1 1 -1 -1 1 -1];

% Thedata descripted in the chapter
initialData = [1 -1 1 1 -1 -1 -1 1 mark 1 -1 1];

%Results of the correlation
correlation = zeros(1, length(data));

% The correlation is performed multiplying all the vectors with all the
%   shifts. The shifts will be done in data
data = initialData;
for x = 1:length(data)
    correlation(x) = mark*data(1:length(mark))'/length(mark);
    data = circshift(data, -1);
end

stem(0:length(correlation)-1, correlation);

% shiftSync = find(