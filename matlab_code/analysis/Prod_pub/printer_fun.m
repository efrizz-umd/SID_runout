% printerfun.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [L] = printer_fun(fileID,print_data,printname)

% ************************************************************************
% This function prints double arrays to text
%   print out a variable in an array style given
%   only works with single column/vector arrays
%   file must already have been opened outside function
%
% % ----------- output ----------- %
% L - length of print data (unused)
% % ----------- intput ----------- %
% - fileID - file to write to
% - print_data - array to print
% - printname - in the text file, will be printname = print_data
%
% ************************************************************************



% build the strings to print
L = length(print_data);
print_string = strcat(printname,' = [');

for jj = 1:L
    if jj < L
        print_string = strcat(print_string, num2str(print_data(jj)),", ");
    else
        print_string = strcat(print_string, num2str(print_data(jj)),"];\n");
    end
end

%print_string = strcat(print_string,"];\n");

fprintf(fileID,print_string);

end
