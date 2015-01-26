
function [data, locs, markers, chrs] = read_data( list_files, info )

%list_files = { 'worm_chrI.txt', 'worm_chrII.txt','worm_chrI.txt' };
%list_spcs = [ 1 1 2  ];

C = length( unique(info.condition ) );

sp = 1;
i_sp = ones(C,1);
data = cell(C,1);
chrs = cell(C,1);
locs = cell(C,1);
for i=1:length(list_files)
	tmp = importdata( list_files{i} );
	sp = info.condition(i);
	data{sp}{i_sp(sp),1} = tmp.data(:,2:end);
	chrs{sp}{i_sp(sp),1} = info.chrs{i};
	locs{sp}{i_sp(sp),1} = tmp.data(:,1);
	i_sp(sp) = i_sp(sp)+1;
	markers = tmp.colheaders;
	if ( isempty( tmp.colheaders{1} ) )
		markers = tmp.colheaders(2:end);
	else
		markers = tmp.colheaders;
	end
 	if (i > 1 &  isequal( pmarkers, markers ) == 0 )
		fprintf( 'marker set different...\n' );
		break;
	end	
	pmarkers = markers;
end

end

