% This is a matlab script with one main function (mutinf) and three subfunctions
% (pairwisemutinf, margent, and jointent).  m-files are named after the main
% function, and subfunctions are only callable from within the m-file.  Besides
% mutinf, pairwisemutinf, margent, and jointent, other functions used in this
% script are core Matlab functions.  You can look up how Matlab functions work using
% the Matlab help system either through Help->MATLAB Help or using the command
% "help functionName".

% When you have completed coding the rest of margent so that it successfully calculates
% the marginal entropy, compute the mutual information of the provided tRNA alignment.
% Use "load" to load the tRNA data into a Matlab matrix.  Then run mutinf (you must
% be in the same directory as the mutinf.m file) on the
% data matrix to return a matrix mutual information.  This data matrix comes from
% a multiple sequence alignment of 1678 tRNA sequences where A, C, G, and U have
% been replaced by 1, 2, 3, and 4, respectively.  This has been done for easy loading
% and processing in Matlab.  Finally, plot the mutual
% information to see which base pairs involve covarying bases.
% The final set of Matlab commands will be:
%   D = load('trna.dat')
%   M = mutinf(D)
%   pcolor(M)
%   colorbar(cool)

% Functions are subroutines in Matlab that receive input parameters and return
% computed values.  Function definitions look like this:
%   function returnValue = functionName(parameter1, parameter2, parameter3)
%     CODE FOR COMPUTING THIS FUNCTION . . .
%   end
% returnValue is the name of the variable returned by the function.  For example,
% in the jointent function, "jentropy" is the returned variable.  When jointent
% is called, it runs through its calculations, and when the final "end" statement
% is reached, the current value of jentropy is returned.  The parameters to
% jointent are the two columns of data needed to compute the joint entropy.

% Programming:
%   length: get the number of elements in an array
%   sum: sum the values in an array or a matrix
%   ~=: relational operator "not equals"; the following code
%          if j ~= 1
%            j = j + 1
%          end
%        will increment j only if j does not equal 1.
% Display commands:
%   pcolor: 2D plot routine
%   colormap(cool): change the color scale used for plotting; makes the
%     underlying grid easier to see
%   surf: 3D plot routine
%   shading interp: interpolate coloring in 3D plot

% Data: data array
function M = mutinf(Data)
  
  [rows,cols] = size(Data);
  M = zeros(cols,cols);
  % Loop through each pair of columns in the alignment
  for i = 1:cols-1
    for j = i+1:cols
      % Calculate the mutual information between two specific columns
      M(i,j) = pairwisemutinf(Data(:,i),Data(:,j));
    end
  end
  
end



% Mutual information between columns i & j
function m = pairwisemutinf(i,j)
  m = margent(i) + margent(j) - jointent(i,j);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPLETE THIS FUNCTION %  it will be very similar to jointent, but less complicated
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marginal entropy (Shannon entropy of column i)
function entropy = margent(i)
  
  % Initialization
  % baseCount is a vector storing counts of A, C, G, & U
  baseCount = [0 0 0 0];
  entropy = 0;

  % Count occurance of each base in this column
  
  % Calculate total number of bases (rows) in the column

  % Loop through the different base types

    % Compute the frequency of a given base (inside loop)

    % Do not calculate log(0); "~=" means "not equal" (inside loop)

      % Add new entropy term to the total entropy (inside loop)

end
  

  
% Joint entropy for columns i & j
function jentropy = jointent(i,j)

  % Initialization
  jentropy = 0;
  basepairCount = [0 0 0 0;
                   0 0 0 0;
                   0 0 0 0;
                   0 0 0 0];

  % Count occurance of each base pair in this pair of columns
  for k = 1:length(i)
    if ( i(k) == 1 )
      if ( j(k) == 1 )
        basepairCount(1,1) = basepairCount(1,1) + 1;
      elseif ( j(k) == 2 )
        basepairCount(1,2) = basepairCount(1,2) + 1;
      elseif ( j(k) == 3 )
        basepairCount(1,3) = basepairCount(1,3) + 1;
      elseif ( j(k) == 4 )
        basepairCount(1,4) = basepairCount(1,4) + 1;
      end
    elseif ( i(k) == 2 )
      if ( j(k) == 1 )
        basepairCount(2,1) = basepairCount(2,1) + 1;
      elseif ( j(k) == 2 )
        basepairCount(2,2) = basepairCount(2,2) + 1;
      elseif ( j(k) == 3 )
        basepairCount(2,3) = basepairCount(2,3) + 1;
      elseif ( j(k) == 4 )
        basepairCount(2,4) = basepairCount(2,4) + 1;
      end
    elseif ( i(k) == 3 )
      if ( j(k) == 1 )
        basepairCount(3,1) = basepairCount(3,1) + 1;
      elseif ( j(k) == 2 )
        basepairCount(3,2) = basepairCount(3,2) + 1;
      elseif ( j(k) == 3 )
        basepairCount(3,3) = basepairCount(3,3) + 1;
      elseif ( j(k) == 4 )
        basepairCount(3,4) = basepairCount(3,4) + 1;
      end
    elseif ( i(k) == 4 )
      if ( j(k) == 1 )
        basepairCount(4,1) = basepairCount(4,1) + 1;
      elseif ( j(k) == 2 )
        basepairCount(4,2) = basepairCount(4,2) + 1;
      elseif ( j(k) == 3 )
        basepairCount(4,3) = basepairCount(4,3) + 1;
      elseif ( j(k) == 4 )
        basepairCount(4,4) = basepairCount(4,4) + 1;
      end
    end
  end

  % Calculate total number of base pairs  
  totalBasepairs = sum(sum(basepairCount));
  [numRows,numCols] = size(basepairCount);
  % Loop through the different base pair types (e.g. A:A, A:C, A:G, etc.)
  for k = 1:numRows
    for m = 1:numCols
      % Compute the frequency of a given base pair
      jointProb = basepairCount(k,m)/totalBasepairs;
      % Do not calculate log(0); "~=" means "not equal"
      if jointProb ~= 0
        % Add new entropy term to the total joint entropy
        jentropy = jentropy - (jointProb * log(jointProb));
      end
    end
  end

end