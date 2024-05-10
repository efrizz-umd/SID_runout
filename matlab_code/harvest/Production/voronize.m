%% voronize.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [vor_x, vor_y, vor_z, vor_vol, vor_r] = voronize(id1,x1,y1,z1,r1,len,height,x_start,x_stop)

% Function takes in position and poutputs the voronoi volume (along with the newly sorted positions and radius)
% calls the voronoi function at the command line
% voro++ must have been compiled and installed
% output of voro++ will be input_filename.vol
% the bounds will need to be adjusted based on the channel dimesions

% Read Positions
delimiter = ' ';

L = length(x1);

% write file to perform voronoi analysis on
% format is file id, x, y, z
format_write = '%d %.12f %.12f %.12f %.12f\r\n';
filename_write = "./voro";
filename_write = convertStringsToChars(filename_write);
filename_write = strcat(filename_write,'me.txt');
fileID_write = fopen(filename_write,'w');

for jj = 1:L
    fprintf(fileID_write,format_write,[id1(jj), x1(jj),y1(jj),z1(jj),r1(jj)]);
end
fclose(fileID_write);

% assign strings - removing slight buffer so that we don''t have to tell voro about walls
l_s = num2str(len - 0.01);
h_s = num2str(height);
f_s = num2str(0.01); % floor
x_m = num2str(x_start);
x_p = num2str(x_stop);

% run voro++ at the command line
% voro++ [opts] <x_min> <x_max> <y_min> <y_max> <z_min> <z_max> <filename>
% -px periodic in x
% the shrunken volume avoids end wall interference and also negates the
command = ['voro++ -r -px -o ' x_m ' ' x_p ' 0.01 ' l_s ' 0.01 ' h_s ' ./vorome.txt'];
system(command);


% now we need to read in the output and get the voronoi volume and
% corresponding coordinates
file_voro = './vorome.txt';
file_voro = convertStringsToChars(file_voro);
file_voro = strcat(file_voro,'.vol');
format_voro = '%f %f %f %f %f%[^\n\r]';

% open the file
fileID_voro = fopen(file_voro,'r');
% read in the newly created file containing each particle''s voronoi volume
dataArray_voro = textscan(fileID_voro, format_voro, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
% close the file
fclose(fileID_voro);

% assign output variables based on read data array
vor_x = dataArray_voro{:, 2};
vor_y = dataArray_voro{:, 3};
vor_z = dataArray_voro{:, 4};
vor_vol = dataArray_voro{:, 5};
vor_r = dataArray_voro{:,6};


end
