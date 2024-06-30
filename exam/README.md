I have implemented the Gauss-Newton method for minimization, which is made to minimize a sum of squares through some vector c.
I have completed three tasks, which I believe serve well as A, B, and C tasks as in the homeworks.

A
Implementing the Gauss-Newton algorithm and applying it to Rosenbrock's valley function and Himmelblau's function, which can both be conveniently written as the sum of two squares.
See Out.txt.

B
For my second task I fitted the PSF (point spread function) of white dwarf star to a 2D Gaussian using my Gauss-Newton method.

2D Gauss:
f(x, y) = A· exp( - ( (x-x0)² / (2σx²) + (y-y0)² / (2σy²) ) )

The data comes from one of my own observations using FUT (det Fjernstyrede UndervisningsTeleskop), which targets a nearby star, LAWD 37, which I am working on measuring the paralax of. The datafile for the observation is a 4096x4096 image, which can be represented as a matrix of values, where a higher value corresponds to more light in that pixel. Using python, I extrated the pixels in a rectange containing my target star, and wrote it to a file (LAWD\_37.fits => LAWD\_37.data).

For the error in each pixel, I used the square root, since it is Poisson distributed (roughly).

I then created separate equations for each pixel, that take the coordinates and a vector with the 2D Gauss function parameters, and give Chi value (r\_ij(c) - data[i, j]) / Sqrt(data[i, j]).

Lastly, I made my Gauss-Newton function also calculate the covariance matrix, so I can get the uncertainties of my fitting parameters.

The results of my fit are in Out.txt. I have also plotted the data and the fit in 3D in LAWD\_37.svg.


C
Lastly, I wanted to time my method with respect to N, when fitting a simple function,
f(x) = c[0] * Sin(c[1] * x + c[3]) * Cos(c[2] * x + c[3]).

I created N points in a defined interval, added noise following a normal destribution, and fitted to these points through the four constants in the function.
For the N = 100 case, I have written the results the fit to Out.txt, and shown the fited function and the data with noise in Simple.svg.

Timing.svg shows the execution time with respect to N. It appears to be linear.
