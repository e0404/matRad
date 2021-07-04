[![View Multivariate Polynomial Regression on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/34918-multivariate-polynomial-regression)

<html><body><div class="content"><h1>Example For MultiPolyRegress</h1><p>X is your Data matrix. 500 data point with 5 dimensions. Another way to look at this is 500 samples of 5 independent variables. Y is your observation vector 500 by 1. You want to find a good polynomial fit of columns of X to Y. Lets say you decided fit a 2nd degree polynomial to all 5 independent variables. And you are for the moment, interested in fitting the standard polynomial basis without further meddling with the terms.</p><h2>Contents</h2><div><ul><li><a href="#1">How to Use the Inputs</a></li><li><a href="#2">Plain</a></li><li><a href="#3">Normalization - Range</a></li><li><a href="#4">Figure</a></li><li><a href="#5">PV</a></li><li><a href="#6">How to Use the Outputs</a></li><li><a href="#7">PowerMatrix</a></li><li><a href="#8">Scores</a></li><li><a href="#9">Polynomial Expression</a></li><li><a href="#10">Cofficients</a></li><li><a href="#11">Legend</a></li><li><a href="#12">yhat</a></li><li><a href="#13">Residuals</a></li><li><a href="#14">Goodness of Fit Measures</a></li></ul></div><h2>How to Use the Inputs<a name="1"></a></h2><h2>Plain<a name="2"></a></h2><pre class="codeinput">load <span class="string">Example.mat</span>
reg=MultiPolyRegress(X,Y,2) <span class="comment">% Gives you your fit.</span>
</pre><pre class="codeoutput">
reg = 

           FitParameters: '-----------------'
             PowerMatrix: [21x5 double]
                  Scores: [500x21 double]
    PolynomialExpression: [21x2 table]
            Coefficients: [21x1 double]
                    yhat: [500x1 double]
               Residuals: [500x1 double]
           GoodnessOfFit: '-----------------'
                 RSquare: 0.9392
                     MAE: 0.0334
                  MAESTD: 0.0481
           Normalization: '1-to-1 (Default)'
      LOOCVGoodnessOfFit: '-----------------'
               CVRSquare: 0.9280
                   CVMAE: 0.0366
                CVMAESTD: 0.0590
         CVNormalization: '1-to-1 (Default)'

</pre><h2>Normalization - Range<a name="3"></a></h2><p>Different error definition ONLY in the calculation of MAE, MAESTD, CVMAE and CVMAESTD. Does not effect the fit.</p><pre class="codeinput">reg=MultiPolyRegress(X,Y,2,<span class="string">'range'</span>)
</pre><pre class="codeoutput">
reg = 

           FitParameters: '-----------------'
             PowerMatrix: [21x5 double]
                  Scores: [500x21 double]
    PolynomialExpression: [21x2 table]
            Coefficients: [21x1 double]
                    yhat: [500x1 double]
               Residuals: [500x1 double]
           GoodnessOfFit: '-----------------'
                 RSquare: 0.9392
                     MAE: 0.0293
                  MAESTD: 0.0313
           Normalization: 'Range'
      LOOCVGoodnessOfFit: '-----------------'
               CVRSquare: 0.9280
                   CVMAE: 0.0313
                CVMAESTD: 0.0346
         CVNormalization: 'Range'

</pre><h2>Figure<a name="4"></a></h2><p>You would like to see a scatter plot of your fit.</p><pre class="codeinput">reg=MultiPolyRegress(X,Y,2,<span class="string">'figure'</span>);
</pre><img vspace="5" hspace="5" src="http://ahmetcecen.github.io/MultiPolyRegressExample/Example_01.png" alt=""> <h2>PV<a name="5"></a></h2><p>You would like to limit the observed powers of certain terms in your polynomial. For example, you do not want the 1st and 4th Independent Variables (x1 and x4) to have second order terms (x1^2 or x4^2). Notice you have to explicitly write how high each term can go in powers, so I would also state I am fine with (x2 x3 and x5) having 2nd order terms.</p><pre class="codeinput">reg=MultiPolyRegress(X,Y,2,[1 2 2 1 2]);
PolynomialFormula=reg.PolynomialExpression
</pre><pre class="codeoutput">
PolynomialFormula = 

    Coefficient     Term  
    ___________    _______

      0.0022642    'x5'   
      0.0058919    'x4'   
    -0.00049119    'x4.x5'
        0.01644    'x3'   
    -0.00098813    'x3.x5'
     7.6129e-05    'x3.x4'
       0.014969    'x2'   
     -0.0023337    'x2.x5'
      0.0028077    'x2.x4'
    -0.00012646    'x2.x3'
      -0.027613    'x1'   
    -0.00036617    'x1.x5'
    -0.00043459    'x1.x4'
    -0.00011518    'x1.x3'
     -0.0009348    'x1.x2'
         3.9964    ''     
      0.0004941    'x2^2' 
     0.00014775    'x3^2' 
       0.010017    'x5^2' 

</pre><h2>How to Use the Outputs<a name="6"></a></h2><pre class="codeinput">reg=MultiPolyRegress(X,Y,2);
</pre><h2>PowerMatrix<a name="7"></a></h2><p>You have a new data point you would like to evaluate using the computed fit. Lets assume for the sake of argument that the 250th row of X is in fact a new data point.</p><p>Unless you have a stake in deeply understanding this code, don't try to make sense of the NewScores matrix, or what follows. I sometimes have to stare at it for a couple minutes to figure it out myself. I am happy discuss this in detail upon specific request.</p><p>You have to repeat this procedure for every new data point. It might be time saving to write a function that does this automatically, however I never needed this functionality, so I wouldn't count on me writing that.</p><pre class="codeinput">NewDataPoint=X(250,:);
NewScores=repmat(NewDataPoint,[length(reg.PowerMatrix) 1]).^reg.PowerMatrix;
EvalScores=ones(length(reg.PowerMatrix),1);
<span class="keyword">for</span> ii=1:size(reg.PowerMatrix,2)
    EvalScores=EvalScores.*NewScores(:,ii);
<span class="keyword">end</span>
yhatNew=reg.Coefficients'*EvalScores <span class="comment">% The estimate for the new data point.</span>
</pre><pre class="codeoutput">
yhatNew =

    5.2877

</pre><h2>Scores<a name="8"></a></h2><p>Unless you have a stake in deeply understanding this code, don't try to make sense of the Scores matrix, chances are you won't ever need to use it.</p><h2>Polynomial Expression<a name="9"></a></h2><p>You would like to see the actual formula of the fit,</p><pre class="codeinput">PolynomialFormula=reg.PolynomialExpression
</pre><pre class="codeoutput">
PolynomialFormula = 

    Coefficient     Term  
    ___________    _______

      0.0052679    'x5'   
      0.0073888    'x4'   
    -8.7941e-05    'x4.x5'
       0.016723    'x3'   
    -0.00097694    'x3.x5'
     8.3902e-05    'x3.x4'
       0.015417    'x2'   
     -0.0025415    'x2.x5'
       0.002392    'x2.x4'
    -0.00018939    'x2.x3'
      -0.028576    'x1'   
    -0.00045571    'x1.x5'
    -0.00037732    'x1.x4'
    -0.00010521    'x1.x3'
    -0.00047654    'x1.x2'
         4.0108    ''     
    -4.9811e-05    'x1^2' 
    -0.00026949    'x2^2' 
     0.00014575    'x3^2' 
    -6.5765e-05    'x4^2' 
       0.011599    'x5^2' 

</pre><h2>Cofficients<a name="10"></a></h2><p>This was shown earlier at the input examples.</p><h2>Legend<a name="11"></a></h2><p>This was shown earlier at the input examples.</p><h2>yhat<a name="12"></a></h2><p>This is the vector of estimates using your new fit. The scatter plot you see when you use the 'figure' option is generated using scatter(yhat,y).</p><h2>Residuals<a name="13"></a></h2><p>This is defined as y-yhat. Can be used for a residual plot to see if ordinary least squares assumptions hold true.</p><h2>Goodness of Fit Measures<a name="14"></a></h2><p>These are useful not only in assesing the accuracy of your fit, but also comparing different candidates. For example, lets see how different powers compare for the same fit for the above dataset. I personally would like to use CVMAE as my comparative error measure, since it is more sensitive to overfitting.</p><p>It turns out, the second degree polynomial is the best option. One way to interpret this number is saying the fit makes in average a 3.66% error with respect to the original Y when estimating.</p><pre class="codeinput"><span class="keyword">for</span> ii=1:5
    reg=MultiPolyRegress(X,Y,ii);
    CVMAE(ii)=reg.CVMAE;
<span class="keyword">end</span>
CVMAE
</pre><pre class="codeoutput">
CVMAE =

    0.0383    0.0366    0.0416    0.1255    3.7904

</pre><p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R2014b</a><br></p></div></body></html>
