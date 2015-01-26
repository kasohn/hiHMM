% This script is used for convert mat file to .txt file.
% The .txt files were used in this study to compare power between hiHIMM and previous methods (ChromHMM, Segway).
input_list = glob('simulation_*_mat/*.mat');

for input = input_list'
	file_name = input{1};
	load(file_name);
	inp_file_name = strrep( file_name, '_mat', '_txt' );
	inp_file_name = strrep( inp_file_name, '.mat', '_signal.txt');

	YY = [ Y{1}{1}' ; Y{2}{1}' ; Y{3}{1}' ];

	fid = fopen( inp_file_name, 'w' );
	fprintf( fid, 'Cell\tchr1\n' );

	for i = 1:8
		if i > 1
			fprintf( fid, '\t');
		end
		fprintf( fid, 'Mark%d', i );
	end
	fprintf( fid, '\n' );
	fclose(fid);

	dlmwrite( inp_file_name, YY, 'delimiter', '\t', '-append' );

	S = [ Strue{1}{1}' ; Strue{2}{1}' ; Strue{3}{1}' ];

	true_file_name = strrep( file_name, '_mat', '_txt' );
	true_file_name = strrep( true_file_name, '.mat', '_true.txt');

	dlmwrite( true_file_name, S );
end
