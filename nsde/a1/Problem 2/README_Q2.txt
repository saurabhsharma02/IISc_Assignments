ASSIGNMENT 1
PROBLEM 2 -SAURABH SHARMA
1. Code:
Wrote separate codes for three methods- Euler's Explicit method (Question 2b1.cpp), Second Order Adam Bashforth (Question 2b2.cpp) and 
Fourth Order Runge Kutta method (Question 2b3.cpp).

2. Input File:
Not created in this example as already step sizes are given.

3. Make File: 
Run the files by using their respective commands- "make t1.o" for Euler's Explicit method, 
"make t2.o" for Second Order Adam Bashforth and 
"make t3.o" for Fourth Order Runge Kutta method OR 
to run all three methods, use "make all".

4. MATLAB FILE: (MATLAB_CODE_Q2)
Contains analytical solution plots and function plots of all three methods for all the given step sizes- 0.5,0.25,0.125,0.0625 and 0.03125.
In MATLAB Code, separate arrays for maximum and mean errors are created for each method.
mEE= Maximum Euler Error
aEE= Average Euler Error
mRK= Maximum Runge Kutta Error
aRK= Average Runge Kutta Error
mAB= Maximum Adam Bashforth Error
aAB= Average Adam Bashforth Error
To compare all three methods, separate plots are made for maximum and average error.