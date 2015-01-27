function write_bedfile( S, loc, chrs, col, model, K, odir, tail, binsize )
% 
% loc: center position of each bin 

  fpa = fopen( sprintf('%s/hihmm.model%d.K%d.%s.bed', odir, model, K, tail), 'w' );
  fprintf( fpa, 'track name=model%dK%d%s description="Chromatin States" useScore=0 itemRGB="On"\n', model, K, tail );
  for ii=1:length(S)
        Sc = S{ii};
        bp = find( Sc(2:end) - Sc(1:end-1) ~= 0 );
        ed = reshape( loc{ii}([bp end])+ceil(binsize/2), 1, length(bp)+1); % block right end coordinate 
        st = [ loc{ii}(1)-ceil(binsize/2)+1 ed(1:end-1)+1]; % block left end coordinate 
	idx_bd = [ 1 bp+1 ];
        for j=1:length(st)
                k = Sc( idx_bd(j) );
                fprintf( fpa, '%s\t%d\t%d\t%d\t1000\t+\t%d\t%d\t%d,%d,%d\n', chrs{ii}, st(j), ed(j), k, st(j), ed(j), col(k,1), col(k,2), col(k,3) );
        end
  end
  fclose(fpa);

end


