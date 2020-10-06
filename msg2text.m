function text = msg2text(msg)
% Converts a bit-stream into an ASCII character string
assert(~rem(length(msg), 8), ...
    'The message doesn''t contain an integer number of bytes'); 
text = char(bin2dec(reshape(num2str(msg), 8, []).').');
end

