function writehmm( pi, mu, p0, fn, odir, markers )

%subid = 1;
%mkdir(odir);
%markers= {'H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K79me2','H3K9me3','H4K20me1'};

K = size( pi, 1 );
fp = fopen( [ odir fn '.transition.csv' ], 'w' );
for kk=1:K-1
	fprintf( fp, '"State %d", ', kk );
end
fprintf( fp, '"State %d"\n', K );
for ee=1:size(pi,3)
	pie = pi(:,1:end-1,ee);
	pie = pie ./ ( sum( pie, 2) * ones(1, K)  );
	pc = p0(ee);
	pie = pc * pie + (1-pc) * eye( size(pie,1) );
	for kk=1:K
		fprintf( fp, '%f, ', pie(kk,1:end-1) );
		fprintf( fp, '%f\n', pie(kk,end) );
	end
end
fclose(fp);

fp = fopen( [ odir fn '.emission.csv' ], 'w' );
for kk=1:length(markers)-1
	fprintf( fp, '%s, ', markers{kk} );
end
fprintf( fp, '%s\n', markers{length(markers)} );
for ee=1:size(mu,3)
	mue = mu(:,:,ee)';
	for kk=1:K
		fprintf( fp, '%f, ', mue(kk,1:end-1) );
		fprintf( fp, '%f\n', mue(kk,end) );
	end
end
fclose(fp);


