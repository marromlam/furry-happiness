function [ U ] = URand( N, varargin )
%function [ U ] = URand( N, a, b, M )
%   Generates random uniform distributed numbers by congruence method. 
%
%      IN:       N - number of random numbes to be generated.
%      OUT:      U - vector of uniformly distributed N numbers.


% Congruence method parameters
  M = 2^(32) -1;                                                  % default
  a = 17000 ; b = 0;                                              % default
  if nargin > 1
      a = varargin{1};
  end
  if nargin > 2
      b = varargin{2};
  end
  if nargin > 3
      M = varargin{3};
  end
  
% First element
  U  = zeros((N+9),1);
  U(1) = round(sqrt(N));                                      % for example

% Loop over N+9 elements, first 10 will not be taken in account.
  for l = 2:1:(N+9)
    U(l) = mod(a*U(l-1)-b,M);
  end
  
% Normalization
  U = U(10:(N+9))/M ;

  return                                                       % Stop here!
  
% Check they are 'really' random
  X = U(1:2:length(U)-1);
  Y = U(2:2:length(U));

  plot(X,Y,'o')
  X = U(1:3:length(U)-1);
  Y = U(2:3:length(U));
  Z = U(3:3:length(U));
  plot3 (X,Y,Z,'o')

end

