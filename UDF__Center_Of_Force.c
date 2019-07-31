/* only for double precision (double variable dclaration) */

# include "udf.h"

#if !RP_NODE
	void PRINT_TO_CONSOLE(double A[])
	{
		/* PRINT RESULTS to console */
		Message("Total:\n %10.7lf %10.7lf %10.7lf\n",
		A[0], A[1], A[2]
		);

		return;
	}
#endif


void INPUT_DATA_function(int INPUT_data[],  double *Fixed_coord_value, int thread_id[], int *length_list_of_threads)
{

	FILE *INPUT_file;

	int Symmery_plane, Fixed_coord;
	int i;


	INPUT_file = fopen("INPUT_file.txt","r");

	fscanf(INPUT_file, "%*[^\n]\n");
	fscanf(INPUT_file, "%d\n", &Symmery_plane);
	fscanf(INPUT_file, "%*[^\n]\n");
	fscanf(INPUT_file, "%d\n", &Fixed_coord);
	fscanf(INPUT_file, "%*[^\n]\n");
	fscanf(INPUT_file, "%lf\n", Fixed_coord_value);
	fscanf(INPUT_file, "%*[^\n]\n");


	for(i = 0; ; ++i) {
		// assign the read value to variable (thread_id[i]) and check condition
		if(fscanf(INPUT_file, "%d\n", &thread_id[i]) == 1)
		{
			(*length_list_of_threads)++;
		}
		else
		{
			// if EOF is returned by fscanf, or in case of error, break loop
			break;
		}
	}

	fclose(INPUT_file);


	/* INPUT_data = [1_coord_of_symmetry_plane ; 2_coord_of_symmetry_plane ; 3_coord ; fixed_coord (on symmetry plane)] */

	switch(Symmery_plane){
		case 1: /*  if  Symmery_plane == 1 'XY' => 'INPUT_data=[0,1,2,..]' =>  z = 0  */
			INPUT_data[0]= 0;
			INPUT_data[1]= 1;
			INPUT_data[2]= 2;
			INPUT_data[3]= INPUT_data[Fixed_coord-1];
			break;
		case 2: /*  if  Symmery_plane == 2 'YZ' => 'INPUT_data=[1,2,0,..]' =>  x = 0  */
			INPUT_data[0]= 1;
			INPUT_data[1]= 2;
			INPUT_data[2]= 0;
			INPUT_data[3]= INPUT_data[Fixed_coord-1];
			break;
		case 3: /*  if  Symmery_plane == 3 'ZX' => 'INPUT_data=[2,0,1,..]' =>  y = 0  */
			INPUT_data[0]= 2;
			INPUT_data[1]= 0;
			INPUT_data[2]= 1;
			INPUT_data[3]= INPUT_data[Fixed_coord-1];
			break;
		default:
			break;
	}

	return;
}



DEFINE_ON_DEMAND(center_of_pressure)
{


Domain *domain = Get_Domain(1);  /* the pointer to the domain (? 1 is deaflt for mixtures??????) */
Thread *t = NULL; /* inicialization of a pointer to the face's thread */


double Force[ND_ND], Force_total[ND_ND]; /* ND_ND - is defined as 2 for 2D case and as 3 for 3D */

double origin_0[ND_ND], Moment[ND_ND], Moment_total[ND_ND], Moment_total_CoP[ND_ND];

int i, n, length_list_of_threads, thread_n_id;
int thread_id[200];

int INPUT_data[4];
double Fixed_coord_value;


double Center_of_pressure[3];


length_list_of_threads = 0;
Fixed_coord_value = 0.0;

/* READ INPUT FILE*/

INPUT_DATA_function(INPUT_data, &Fixed_coord_value, thread_id, &length_list_of_threads);

/* START CALCULATE CENTER OF PRESSURE */

for(i=0 ; i<ND_ND ; i++)
{
	origin_0[i]=0.0;
	Center_of_pressure[i] = 0.0;

	Force[i] = 0.0;
	Force_total[i]  = 0.0;

	Moment[i] = 0.0;
	Moment_total[i] = 0.0;
	Moment_total_CoP[i] = 0.0;
}


for(n=0 ; n<length_list_of_threads ; n++)
{
	thread_n_id = thread_id[n];
	t = Lookup_Thread(domain,thread_n_id);

	Compute_Force_And_Moment(domain, t, origin_0, Force, Moment, TRUE);

	  /* This function needs size of 3 for force and moment.  It works in parallel and
	the common arguments are the of the same type as your arguments.
	In addition you have to pass the domain as first arg and a boolean as
	last argument. This boolean has to be TRUE if you also call the function
	on the host else it has to be FALSE.  */

	for(i=0 ; i<ND_ND ; i++)
	{
		Force_total[i] += Force[i];
		Moment_total[i] += Moment[i];
	}
}

Center_of_pressure[INPUT_data[3]] = Fixed_coord_value;

if(INPUT_data[3] == 1) /* Fixed_coord ; fixed  INPUT_data[0] */
{
	Center_of_pressure[INPUT_data[1]] = ( Force_total[INPUT_data[1]] / Force_total[INPUT_data[0]] ) * Fixed_coord_value - Moment_total[INPUT_data[2]] / Force_total[INPUT_data[0]];
	/* DODAC ERROR JAK DZIELENIE PRZEZ 0*/
}
else if(INPUT_data[3] == 2) /* Fixed_coord ; fixed  INPUT_data[1] */
{
	Center_of_pressure[INPUT_data[0]] = ( Force_total[INPUT_data[0]] / Force_total[INPUT_data[1]] ) * Fixed_coord_value + Moment_total[INPUT_data[2]] / Force_total[INPUT_data[1]];
}

/* END CALCULATE CENTER OF PRESSURE */


for(n=0 ; n<length_list_of_threads ; n++)
{
	thread_n_id = thread_id[n];
	t = Lookup_Thread(domain,thread_n_id);
	Compute_Force_And_Moment(domain, t, Center_of_pressure, Force, Moment, TRUE);

	for(i=0 ; i<ND_ND ; i++)
	{
		Moment_total_CoP[i] += Moment[i];
	}
}

/* PRINT RESULTS to CONSOLE */

#if !RP_NODE

	Message("FORCES:\n");
	PRINT_TO_CONSOLE(Force_total);

	Message("MOMENTS(0,0,0):\n");
	PRINT_TO_CONSOLE(Moment_total);


	Message("Coords of center of pressure: [%20.17lf,%20.17lf,%20.17lf]\n",
	Center_of_pressure[0], Center_of_pressure[1], Center_of_pressure[2]);

	Message("Moments for center of pressure: [%10.7lf,%10.7lf,%10.7lf]\n",
	Moment_total_CoP[0], Moment_total_CoP[1], Moment_total_CoP[2]);

#endif

}

/********** END OF UDF *****************/