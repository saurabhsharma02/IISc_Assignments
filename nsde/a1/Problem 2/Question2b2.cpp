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
    long double k1, k2, k;
    long double sze[] = {0.5,0.25,0.125,0.0625,0.03125};
    for (long double h : sze)
    {

        //Defining and declaring number of nodes which is integer value where L is number of nodes.

        int L = static_cast<int>(2/h) + 1; 
        std::vector<long double> t(L);

        for (int i = 0 ; i < L ; i++)
            t[i] = i*h;

        std::vector<long double> y(L);

         y[0] = 1;
        //Runge Kutta Second Order Method to caculate first element
        k1= h*f (t[0],y[0]);
        k2= h*f (t[0]+h,y[0]+k1);
        k= k1+k2;
        y[1]= y[0]+k;

        
        //Adam Bashforth Seond order Method "for" loop.
        for (int i = 2 ; i < L ; i++)
        {
            y[i] = y[i-1] + h*f(t[i-1],y[i-1]);
        }


        //fstream to create files and write data obtained from Second Order Adam Bashforth Method.
        
        std::fstream f;
        f.open("AB2b2_"+std::to_string(h)+".csv",std::ios::out);

        for(int i = 0 ; i < L ; i++)
            f << t[i] << "," << y[i] << "\n";
        f.close();
    }

}
