#include<fstream>
#include<vector>
#include<cmath>


//Defining function as the differential equation given in the problem.

//f= dy/dt = y*t^2 -1.1*y

long double f (long double t, long double y)
{
    long double dydt= y*t*t -1.1*y;
    return dydt;
}

int main()
{
    long double sze[] = {0.5,0.25,0.125,0.0625,0.03125};
    for (long double h : sze)
    {

        //Defining and declaring number of nodes which is integer value where L is number of nodes.

        int L = static_cast<int>(2/h) + 1; 
        std::vector<long double> t(L);

        for (int i = 0 ; i < L ; i++)
            t[i] = i*h;

        std::vector<long double> y(L);

        //Explicit Euler Method "for" loop.

        y[0] = 1;
        for (int i = 1 ; i < L ; i++)
        {
            y[i] = y[i-1] + h*f(t[i-1],y[i-1]);
        }


        //fstream to create files and write data obtained from Euler Explicit Method.
        
        std::fstream f;
        f.open("Euler2b1_"+std::to_string(h)+".csv",std::ios::out);

        for(int i = 0 ; i < L ; i++)
            f << t[i] << "," << y[i] << "\n";
        f.close();
    }

}
