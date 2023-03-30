function [ pf, varargout ] = NonLinearModelFit( model, x , y, varargin )
%function [ pf, varargout ] = NonLinearModelFit( model, x , y, varargin )
%  Performs a non-linear fit of x and y data to given model(x,p) by
%  minimizing chi2-sum-of-residuals over p_i parameters of the model (only
%  up to 10 parameters are allowed). It provides confidence intervals using
%  Fisher tests and ANOVA analysis and stores fit over data plot in picture
%  FitFigure.eps. Only f, x and y are compulsury, other arguments are
%  optional and must be provided in the following order:
%
%      IN:   model - function model(x,p) with .*, ./ etc operators which
%                    depends on x (X axis data) and p(i) parameters.
%               p0 - seek to which fit parameters must be near, it is not 
%                    compulsury but if model is highly non-linear, fit may 
%                    not converge otherwise (default is 1s).
%               cl - confidence level of parameters confidence interval, it
%                    must be between 0 and 1 (default is 0.95).
%        fitfigure - title and labels struct varaible; it must have writen
%                    in LaTeX notation:
%                    fitfigure.titlestr = 'title' (default is empty),
%                    fitfigure.xlabelstr = 'x label' (default is x(u)),      
%                    fitfigure.ylabelstr = 'y label' (default is y(u)).
%          realtol - tolerance to take decision if parameters are 'really'
%                    complex, it is about the natural number that defines
%                    the 1e(realtol) order of magnitude (default is 2).
%      OUT:     pf - final parameters.
%             Pars - interval of cl confidence level of the parameters, it 
%                    is an out second argument of the function.


% Decreasing model inputs and backup original model
  plotmodel = @(x,p) model(x,p);
  model     = @(p)   model(x,p);

% Get number of parameters (NoP) of an inline function.
  check = 0; NoP = 1; 
  while (check == 0) && (NoP < 10)
    try
      model(ones(NoP,1));
      check = 1;
    catch
      check = 0; NoP = NoP + 1;
    end
    if NoP == 10
        error(['Check your model parameters and try again. Remember',...
               'this function can not fit a >10 parameter model.'])
    end
  end

% Set parameters  
  p0 = ones(NoP,1);                             % First guess (default one)
  f = @(p) sum((y - model(p)).^2);                          % Chi2 function
  realtol = 2;               % 1e(realtol), tolerance about complex numbers
  cl = 0.95;                               % Confidence level (default one)
  fitfigure.titlestr = '';                      % FitFigure title (default)
  fitfigure.xlabelstr = '$x\, (\mathrm{u})$'; % FitFigure X label (default)
  fitfigure.ylabelstr = '$y\, (\mathrm{u})$'; % FitFigure Y label (default)
  if nargin > 3
    p0 = varargin{1}; 
    if length(p0) ~= NoP
        error('Seek size and number of parameters are inconsistent.')
    end
  end
  if nargin > 4
    cl = varargin{2};
    if (cl > 1) || (cl < 0)
        error(['Confidence level is a proability, and it must be ',...
               'between 0 and 1.'])
    end
  end
  if nargin > 5 
    try
      fitfigure = varargin{3};
    catch
      warning('Figure labels are not properly setted. Using default ones.')
    end
  end
  if nargin > 6
    realtol = varargin{4}; 
  end
  h = min(sqrt(eps)*p0);              % Chi2 numerical differentiation step 

% Creating df1 and df2 to calculate df derivative at each p_i (aka w_i)
  for l = 1:1:NoP
    str1 = ''; str2 = '';
    for m = 1:1:l-1
      if m == 1
        str1 = 'w('+string(m)+'),';
      else
        str1 = str1 + 'w('+string(m)+'),';
      end
    end
    for m = l+1:1:NoP
      if m == l+1
        str2 = ',w('+string(m)+')';
      else
        str2 = str2 + ',w('+string(m)+')';
      end
    end
    if l == 1
      sf1 = strcat('f([',str1,'w(',string(l),')+h',str2,'])');
      sf2 = strcat('f([',str1,'w(',string(l),')-h',str2,'])');
    else 
      sf1 = strcat(sf1,';','f([',str1,'w(',string(l),')+h',str2,'])');
      sf2 = strcat(sf2 ,';','f([',str1,'w(',string(l),')-h',str2,'])');
    end
  end
  eval(strcat('df1 = @(w)[',sf1,'];'))              % f_i vector at p_i + h
  eval(strcat('df2 = @(w)[',sf2,'];'))              % f_i vector at p_i - h

% Using Newton Raphson to find zeros of system of NoP equations.
  [pf] = NewtonRaphson(@(w) (df1(w) - df2(w))/(2*h), p0, 5e2);
  if ( (NumPow(imag(pf))) <= (NumPow(real(pf))-realtol) )
      [pf] = real(pf);
  else
    error(['Complex value computed by model function, ',...
           'fitting cannot continue. ']);
  end

% Plot model and data.
  set(0,'defaultTextInterpreter','latex')          % LaTeX Text Interpreter
  set(gca,'TickLabelInterpreter', 'latex');
  fig = figure(1); hold all; box on; pbaspect([(1+sqrt(5))/2, 1, 1])
  a = (0.8)^( sign(min(x)))*sign(min(x))*abs(min(x));
  b = (0.8)^(-sign(max(x)))*sign(max(x))*abs(max(x));
  plotx = real(linspace(a,b,1000)');           % Continuous x plot variable
  plot(plotx,real(plotmodel(plotx,pf)),'-b','LineWidth',1.25)  % Model plot
  scatter(x,y,'g','filled')                                     % Data plot
  title(fitfigure.titlestr,'Interpreter','latex')
  xlabel(fitfigure.xlabelstr,'Interpreter','latex')
  ylabel(fitfigure.ylabelstr,'Interpreter','latex')
  print(fig,'FitFigure','-depsc','-tiff')
  
% Calculate parameter errors
  SSopt = real(f(pf)); SStot = sum((y - mean(y)).^2); 
  dof = length(x) - NoP; FisherStat = finv(cl,NoP,dof);
  SS2cl = SSopt*(FisherStat*(NoP/(dof))+1);
  R2 = 1 - SSopt/SStot;
  AR2 = 1 - (1-R2)*(length(x)-1)/(dof-1);
  X2 = chi2inv(cl,dof);
  X2dof = X2/dof;
  [ParMin] = NewtonRaphson( @(w) f(w) - SS2cl, pf-0.1*pf);
  [ParMax] = NewtonRaphson( @(w) f(w) - SS2cl, pf+0.1*pf);
  Pars = sort([ParMin, ParMax],2);
  ParMin = Pars(:,1); ParMax = Pars(:,2);
  
% Table of estimated parameters confidence interval
  rownames = cell(1,NoP);
  for l = 1:1:NoP
    rownames{l} = strcat('p(',string(l),')');
  end
  rownames = cellstr(rownames)';
  colnames = {'Parameter_Minimum','Parameter_Maximun'};
  upf = real(max(pf-ParMin,ParMax-pf));
  if ((NumPow(upf))>=(NumPow(imag(min(ParMax,ParMin)))-realtol))
      ParMin = real(ParMin); ParMax = real(ParMax);
  else
    error(['Complex value computed by model function, ',...
           'fitting cannot continue. Check the results.']);
  end
  
% Print info about the fit
  fprintf('<strong>PARAMETER ESTIMATES</strong>\n')
  fprintf('    <strong>p(%i)</strong>\t\t%.10f\n',([[1:1:NoP]',pf])')
  fprintf('\n<strong>GOODNESS OF FIT</strong>\n')
  fprintf('    <strong>SSE:</strong>\t\t%.10f\n',SSopt)
  fprintf('    <strong>R-square:</strong>\t\t%.10f\n',R2)
  fprintf('    <strong>Adjusted R-square:</strong>\t%.10f\n',AR2)
  fprintf('    <strong>X2 percentile:</strong>\t%.10f\n',X2)
  fprintf('    <strong>X2/dof:</strong>\t\t%.10f\n',X2dof)
  fprintf('    <strong>dof:</strong>\t\t%i\n',dof)
  fprintf('    <strong>Confidence level:</strong>\t%.10f\n',cl)
  table(ParMin,ParMax,'RowNames',rownames,'VariableNames',colnames)
  
% Output variables  
  if nargout>1
    varargout{1} = Pars;
  end
  
end


%% Decimal power of given number

function [ Y ] = NumPow( X )
% Calculate the decimal power of a given number.

  Y = round(log10(abs(X+eps)));
  Y = Y - ((10.^Y)>abs(X+eps));
 
end
