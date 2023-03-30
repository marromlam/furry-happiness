function [ root, varargout ] = NewtonRaphson( f, varargin )
%function [ root, varargout ] = NewtonRaphson( f, varargin )
%  Solves the equation f(x) = 0 by the Newton Raphson method using a
%  relaxed condition: x(k+1) = x(k) - J(f(x(k))\f(x(k)). This makes
%  possible not to calculate the inverse of J, wich is a computational
%  expensive process. Only f is compulsury, other arguments are optional
%  and must be provided in the following order:
%
%      IN:       f - function f(x) of Class 1, to which we will find 
%                    the root.
%               x0 - seek to which the root must be near (default is 1s).
%            maxit - max. number of iterations (default is 1e2).
%            eepsi - epsilon error tolerance (default is 1e-14).
%            edelt - delta error tolerance (default is 1e-14).
%      OUT:   root - root of f(x) = 0.
%                l - number of iteration to achieve given tolerance, it
%                    will allways be less or equal than maxit, if it is
%                    equal tolerance wont be achieved (not by default).
%            froot - f(root), value of f at root, it must be close to zero
%                    if the algorithm performed well (not by default). 


% Get dimesion (NoE) of f(x) inline function.
  check = 0; NoE = 1; 
  while (check == 0)
    try
      f(ones(NoE,1));
      check = 1;
    catch
      check = 0; NoE = NoE + 1;
    end
  end
  
% Set parameter values
  x0 = zeros(NoE,1);                                          % First guess
  maxit = 1e2;                           % Default max number of iterations
  eepsi = 1e-14; edelt = 1e-14;                         % Default tolerance 
  if nargin > 1
    x0 = varargin{1}; 
    if length(x0) ~= NoE
      error('Seek size and number of parameters are inconsistent.')
    end
  end
  if nargin > 2
    maxit = varargin{2};
  end
  if nargin > 3
    eepsi = varargin{3};
  end
  if nargin > 4
    edelt = varargin{4};
  end


% Finding root
  for l = 1:1:maxit
    Jf0 = NumJacobian(f,x0);                          % Jacobian of f at x0
    root = Jf0\(-f(x0)) + x0;      % RELAXED, non relaxed was x0-J^-1*f(x0)
    froot=f(root);
    if norm(root-x0,2) < eepsi || norm(froot,2) < edelt
      break;
    end
    x0 = root; 
  end
  if l == maxit
    warning('Max number of iterations.');
  end


% Output variables.
  if nargout>1
    varargout{1} = l;
  end
  if nargout>2
    varargout{2} = froot;
  end
  
end


%% Jacobian module

function [ J ] = NumJacobian( f, x0 )

  [r,c] = size(f(zeros(length(x0))));
  c = max(r,c); r = max(r,length(x0));
  J = zeros(r,c); 
  for l = 1:1:r
    if x0(l) ~= 0
      h = sqrt(eps)*x0(l);
    else 
      h = 1e-14;
    end
    xh1 = x0; xh1(l) = x0(l) + h; 
    xh2 = x0; xh2(l) = x0(l) - h;
    J(l,:)  = (feval(f,xh1) - feval(f,xh2))./(2*h);
  end
  J = transpose(J);
  
end
