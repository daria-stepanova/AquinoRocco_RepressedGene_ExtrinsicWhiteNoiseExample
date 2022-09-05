#include "Declarations.h"

using namespace std;

double f_func(double t, double x, double beta, double k)
{
	return beta - k*x;
}

double g_func(double t, double x, double D)
{
	return -sqrt(D)*x;
}

int main()
{
	// simple SDE with white noise in the degradation rate
	srand(time(0));
	
	ofstream sde_solution_output;
	sde_solution_output.open("Outputs/output_D1e-5.txt");

	double g = 0.1, k = 5e-4, rho = 10.0, R = 100.0;
	double D = 1e-5; // noise intensity

	double t0 = 0.0, dt = 1e-2, tf = 10000000;
	double tlen = tf/dt;
	double print_counter = tlen/(double)20000.0;

	double beta = g/(1+rho*R);

	double x_mode_theory = beta/(k+D);

	cout<<"Theoretical prediction of the stationary mode is "<<x_mode_theory<<endl;

	double x_new = 0.0, x_prev = 0.0;

	default_random_engine generator;
    normal_distribution<double> st_norm_distribution(0.0,1.0);

	double t_i = 0.0, dW = 0.0, gnbar = 0.0;
	double percent_counter = 0.0;
	double count_p = 0.0;
	for(double i = 0.0; i<tlen-1.0; i++)
	{
		dW = st_norm_distribution(generator)*sqrt(2.0*dt);
		t_i = t0+i*dt;
		gnbar = g_func(t_i,x_prev+g_func(t_i,x_prev,D)*dW,D);
		x_new = x_prev + f_func(t_i,x_prev,beta,k)*dt + 0.5*(g_func(t_i,x_prev,D)+gnbar)*dW;
		x_prev = x_new;
		if(i==count_p*print_counter)
		{
			sde_solution_output<<x_new<<endl;
			count_p++;
		}
		
		if(i/(double)tlen >= percent_counter)
		{
			cout<<percent_counter<<" is done"<<endl;
			percent_counter+=0.1;
			sde_solution_output.flush();
		}
	}

	cout<<"done"<<endl;

	sde_solution_output.flush();
	sde_solution_output.close();


	return 0;
}
