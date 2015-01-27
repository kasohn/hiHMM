function [G, p0s ] = SampleAuxVariable( S, p0s, Pi );

C = length(S);

for cc=1:C 
 for j=1:length(S{cc})
  G{cc}{j} = zeros( 1, length(S{cc}{j}) );
  p0 = p0s(cc);
  for tt=2:length(S{cc}{j})
	prob1 = p0 * Pi( S{cc}{j}(tt-1), S{cc}{j}(tt), cc );	
	prob2 = (1-p0) * ( S{cc}{j}(tt-1) == S{cc}{j}(tt) )  ;
	prob1 = prob1 / (prob1+prob2);
	if ( rand() < prob1 )
		G{cc}{j}(tt) = 1;
	else
		G{cc}{j}(tt) = 0;
	end
  end
 end
 gm = cell2mat(G{cc});
 p0s(cc) = mean( gm(:) );
end


