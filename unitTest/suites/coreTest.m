% =========================================================================
% *** FUNCTION coreTest
% ***
% *** matRad coreTest functions
% ***
% =========================================================================
function [status] = coreTest(k)

  % assign the functions to test
  testfunction_handles = {                        ...
                           @check_some_stuff    
                         };


  numFunctions = length( testfunction_handles );

  if (k<=0)
      status = testfunction_handles;
      return;  % This is used for querying numFunctions.

  elseif (k<=numFunctions)
      status = testfunction_handles{k}();
      status.function = func2str(testfunction_handles{k});

  else
      error('testfunctions:outOfBounds', ...
            'Out of bounds (number of testfunctions=%d)', numFunctions);
  end

end

function [stat] = check_some_stuff()
  stat.description = 'check_some_stuff is running..';
  stat.unreliable = isOctave || isMATLAB(); %FIXME: `width` is inconsistent, see #552

  m = [0 1 1.5 1 -1];
  plot(m,'*-'); hold on;
  plot(m(end:-1:1)-0.5,'x--');

  title({'multline','title'});
  legend({sprintf('multi-line legends\ndo work 2^2=4'), ...
        sprintf('second\nplot')});
  xlabel(sprintf('one\ntwo\nthree'));
  ylabel({'one','° ∞', 'three'});

  set(gca,'YTick', []);
  set(gca,'XTickLabel',{});

end