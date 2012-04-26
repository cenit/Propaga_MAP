/*********************************************************************************
Copyright 2010, 2011, 2012 Stefano Sinigardi, Graziano Servizi, Giorgio Turchetti
*********************************************************************************/


#include "Propaga_MAP.h"



class Elemento_magnetico
{
protected:
	double lunghezza, parametro1, posizione_iniziale;

public:
	void setvalues(double l, double k1, double l_start)
	{
		lunghezza=l,
		parametro1=k1;
		posizione_iniziale=l_start;
	}

	double getvalues(int i)
	{
		if      (i==0)	return lunghezza;
		else if (i==1)	return parametro1; 
		else if (i==2)  return posizione_iniziale;
		else
		{
						return -123456;
		}
	}
};


struct Drift : public Elemento_magnetico
{
	Drift(){}
	Drift(double l)
	{
		lunghezza = l;
	}
//	mappa_(double * xoriginale,  double * xevoluto, double zeta, double * parametriElemento);
	static void mappa_(double *,  double *, double, double *);
};


struct Solenoid : public Elemento_magnetico
{
	Solenoid(){}
	Solenoid(double l, double k1)
	{
		lunghezza = l, 
		parametro1 = k1;
	}
	static void mappa_(double *, double *, double, double *);
};


struct Focusing : public Elemento_magnetico
{
	Focusing(){}
	Focusing(double l, double k1)
	{
		lunghezza = l, 
		parametro1 = k1;
	}
	static void mappa_(double *,  double *, double, double *);
};


struct Defocusing : public Elemento_magnetico
{
	Defocusing(){}
	Defocusing(double l, double k1)
	{
		lunghezza = l, 
		parametro1 = k1;
	}
	static void mappa_(double *,  double *, double, double *);
};


void  Drift :: mappa_(double *Xold,  double *Xnew, double zeta, double *Param)
{
	double D[4][4];
	memset((void *)Xnew, 0, 4*sizeof(double));

	D[0][0] = D[1][1] = D[2][2] = D[3][3] = 1.0;
	D[0][1] = D[2][3] = zeta;
	D[0][2] = D[0][3] = D[1][0] = D[1][2] = D[1][3] = D[2][0] = D[2][1] = D[3][0] = D[3][1] = D[3][2] = 0.0;

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			Xnew[i] += D[i][j] * Xold[j];
}


void  Solenoid :: mappa_(double *Xold,  double *Xnew, double zeta, double *Param)
{
	// Xold[] = {x, x', y, y', z, pz, beta}
	// Param[] = {tipo, posizione_iniziale, lunghezza, campo}

	double beta = Xold[6];
	double sk = sqrt((Param[3]/beta)*(Param[3]/beta));
	double alpha = (Param[2]) * sk;
	double phi = (zeta - Param[1]) * sk;

	double *ret = new double[4];
	double M[4][4], R[4][4];
	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			M[i][j] = R[i][j] = 0.0;

	if (zeta < (Param[1]+Param[2]))
	{
		M[0][0] = R[0][0] = M[1][1] = R[1][1] = M[2][2] = R[2][2] = M[3][3] = R[3][3] = cos(phi);

		M[0][1] = M[2][3] = sin(phi) / sk;
		M[1][0] = M[3][2] = -sin(phi) * sk;
		R[0][2] = R[1][3] = sin(phi);
		R[2][0] = R[3][1] = -sin(phi);
		R[1][0] = -sin(phi) * sk;
		R[1][2] = cos(phi) * sk;
		R[3][0] = -cos(phi) * sk;
		R[3][2] = -sin(phi) * sk;
	}
	else
	{
		M[0][0] = R[0][0] = M[1][1] = R[1][1] = M[2][2] = R[2][2] = M[3][3] = R[3][3] = cos(alpha);

		R[0][2] = R[1][3] = sin(alpha);
		R[2][0] = R[3][1] = -sin(alpha);

		M[0][1] = M[2][3] = sin(alpha) / sk;
		M[1][0] = M[3][2] = -sin(alpha) * sk;
	}

	memset((void *)ret, 0, 4*sizeof(double));

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			ret[i] += M[i][j] * Xold[j];

	memset((void *)Xnew, 0, 4*sizeof(double));

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			Xnew[i] += R[i][j] * ret[j];
}


void  Focusing :: mappa_(double *Xold,  double *Xnew, double zeta, double *Param)
{
	// Xold[] = {x, x', y, y', z, pz, beta}
	// Param[] = {tipo, posizione_iniziale, lunghezza, campo}
	double kF = Param[3];

	double zetaF = zeta - Param[1];

	double sqrtkF = sqrt(kF);
	double skF = sin(sqrtkF * zetaF);
	double ckF = cos(sqrtkF * zetaF);
	double shkF = sinh(sqrtkF * zetaF);
	double chkF = cosh(sqrtkF * zetaF);

	double FOC[4][4];
	memset((void *)Xnew, 0, 4*sizeof(double));

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			FOC[i][j] = 0.0;

	FOC[0][0] = FOC[1][1] = ckF;
	FOC[2][2] = FOC[3][3] = chkF;
	FOC[0][1] = skF / sqrtkF;
	FOC[1][0] = -sqrtkF * skF;
	FOC[2][3] = shkF / sqrtkF;
	FOC[3][2] = sqrtkF * shkF;

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			Xnew[i] += FOC[i][j] * Xold[j];
}


void  Defocusing :: mappa_(double *Xold,  double *Xnew, double zeta, double *Param)
{
	// Xold[] = {x, x', y, y', z, pz, beta}
	// Param[] = {tipo, posizione_iniziale, lunghezza, campo}
	double kD = Param[3];

	double zetaD = zeta - Param[1];

	double sqrtkD = sqrt(kD);
	double skD = sin(sqrtkD * zetaD);
	double ckD = cos(sqrtkD * zetaD);
	double shkD = sinh(sqrtkD * zetaD);
	double chkD = cosh(sqrtkD * zetaD);

	double DEF[4][4];
	memset((void *)Xnew, 0, 4*sizeof(double));

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			DEF[i][j] = 0.0;

	DEF[2][2] = DEF[3][3] = ckD;
	DEF[0][0] = DEF[1][1] = chkD;
	DEF[2][3] = skD / sqrtkD;
	DEF[3][2] = -sqrtkD * skD;
	DEF[0][1] = shkD / sqrtkD;
	DEF[1][0] = sqrtkD * shkD;

	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			Xnew[i] += DEF[i][j] * Xold[j];
}


void create_gnuplot_file(std::string gnuplot_filename, std::string run_name, double *param)
{
	std::ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());

	gnuplot_file << "#!/gnuplot" << std::endl;
	gnuplot_file << "FILE1=\"out_" << run_name << ".ppg\"" << std::endl;
	gnuplot_file << "set terminal postscript eps enhanced colour solid rounded \"Helvetica\" 25" << std::endl;
	gnuplot_file << "set output \"graph_" << run_name << ".eps\"" << std::endl;
	gnuplot_file << "set xlabel \"z (cm)\"" << std::endl;
	gnuplot_file << "plot FILE1 u 4:1 w lines lt 1 lc rgb \"red\" lw 3 t \"x\",\\" << std::endl;
	gnuplot_file << "FILE1 u 4:3 w lines lt 1 lc rgb \"blue\" lw 3 t \"y\"" << std::endl;

	gnuplot_file.close();
}


int main(int argc, char *argv[])
{
	int major_version=1;
	int minor_version=0;
	int fix_release=0;
	std::string release_date="April 25, 2012";
	std::string latest_commit="first stable release";

	std::string run_name="test";
	std::ifstream parametri_sim, parametri_lattice;
	std::ofstream log_file, output_file;
	std::ostringstream log_temp, log_filename, output_filename, gnuplot_filename;
	std::string filename, gnuplot_filename_string, log_temporaneo, parametri_filename;

	bool fallita_lettura_parametri_sim = true, fallita_lettura_parametri_lattice = true;


	for (int i = 1; i < argc; i++)	// * We will iterate over argv[] to get the parameters stored inside.
	{								// * Note that we're starting on 1 because we don't need to know the
									// * path of the program, which is stored in argv[0]
		if (std::string(argv[i]) == "-lattice")
		{
			parametri_lattice.open(argv[i+1]);
			fallita_lettura_parametri_lattice = parametri_lattice.fail();
			log_temp << "Lattice file: " << std::string(argv[i+1]) << std::endl;
			i++;							// so that we skip in the for cycle the parsing of the <nome-file-particelle>
		}
		else if (std::string(argv[i]) == "-params") 
		{
			parametri_sim.open(argv[i+1]);
			fallita_lettura_parametri_sim = parametri_sim.fail();
			parametri_filename = std::string(argv[i+1]);
			log_temp << "Parameters for the simulation in file: " << std::string(argv[i+1]) << std::endl;
			i++;							// so that we skip in the for cycle the parsing of the <nome-file-lattice>
		}
		else if (std::string(argv[i]) == "-out")
		{
			run_name = std::string(argv[i+1]);
			log_temp << "Run name: " << run_name << std::endl;
			i++;
		}
		else
		{
			log_temp << "Invalid argument: " << argv[i] << std::endl;
			return 243;
		}
	}

	log_filename << "LOG_" << run_name << ".ppg";
	log_temporaneo = log_temp.str();
	filename = log_filename.str();
	log_file.open(filename.c_str());
	log_file << "\nPropagaMAP v" << major_version << "." << minor_version << "." << fix_release << "\nRelease date: " << release_date << "\nLatest change: " << latest_commit << std::endl;
	log_file << log_temporaneo;

 	if((fallita_lettura_parametri_sim || fallita_lettura_parametri_lattice ))
	{
		std::cerr << "Missing some required input files: please use -lattice latticefilename" << std::endl;
		std::cerr << "and -params paramsfilename to tell the program the input data, -out name to give" << std::endl;
		std::cerr << "a name to the simulation" << std::endl;
		return 253;
	}

	output_filename << "out_" << run_name << ".ppg";
	filename = output_filename.str();
	output_file.open(filename.c_str());
	gnuplot_filename << "plot_" << run_name << ".plt";
	gnuplot_filename_string = gnuplot_filename.str();


	std::string utile_per_contare, utile_per_leggere;
	int conta_righe_parametri = 0;
	do
	{
		parametri_sim >> utile_per_contare;
		if(parametri_sim.eof()) break;
		parametri_sim.ignore(1000, '\n');
		conta_righe_parametri++;
	}
	while(!parametri_sim.eof());
	parametri_sim.clear();
	parametri_sim.seekg(0,std::ios::beg);

	if (conta_righe_parametri < NUMERO_PARAMETRI_MINIMO_SIMULAZIONE)
	{
		std::cerr << "Missing some required parameter in " << parametri_filename << std::endl;
		return 243;
	}

	double zmax, x0, y0, z0, px0, py0, E0;
	int nsteps;
	const char * parametro_letto;

	for (int i = 0; i < conta_righe_parametri; i++)
	{
		parametri_sim >> utile_per_contare;
		parametri_sim >> utile_per_leggere;
		parametro_letto = utile_per_leggere.c_str();
		if (utile_per_contare == "ZMAX_CM")
		{
			zmax = atof(parametro_letto);
			log_temp << "ZMAX (in cm): " << zmax << std::endl;
		}
		else if (utile_per_contare == "NSTEPS") 
		{
			nsteps = atoi(parametro_letto);
			log_temp << "Number of steps: " << nsteps << std::endl;
		}
		else if (utile_per_contare == "X0") 
		{
			x0 = atof(parametro_letto);
			log_temp << "initial x position (in cm): " << x0 << std::endl;
		}
		else if (utile_per_contare == "Y0") 
		{
			y0 = atof(parametro_letto);
			log_temp << "initial y position (in cm): " << y0 << std::endl;
		}
		else if (utile_per_contare == "Z0") 
		{
			z0 = atof(parametro_letto);
			log_temp << "initial z position (in cm): " << z0 << std::endl;
		}
		else if (utile_per_contare == "PX0") 
		{
			px0 = atof(parametro_letto);
			log_temp << "initial px (normalized): " << px0 << std::endl;
		}
		else if (utile_per_contare == "PY0") 
		{
			py0 = atof(parametro_letto);
			log_temp << "initial py (normalized): " << py0 << std::endl;
		}
		else if (utile_per_contare == "E0") 
		{
			E0 = atof(parametro_letto);
			log_temp << "Energy (in MeV): " << E0 << std::endl;
		}
		else log_temp << "Unrecognized parameter: " << utile_per_contare << std::endl;
	}


	int *N_tipi = new int[N_TIPI_MAGNETICI];
	memset((void *)N_tipi, 0, N_TIPI_MAGNETICI * sizeof(int));
	std::string parametro_1, parametro_2;
	const char * parametro_1s;
	const char * parametro_2s;
	double parametro_letto_1, parametro_letto_2;

	while(1)  // loop infinito
	{
		parametri_lattice >> utile_per_contare;
		parametri_lattice >> parametro_1;
		parametri_lattice >> parametro_2;

		if(parametri_lattice.eof()) break;     // uscita dal loop infinito alla fine del file
		parametri_lattice.ignore(1000, '\n');  // finche' il file non e' finito saltare a fine riga
		if (utile_per_contare == "O") N_tipi[_DRIFT_]++;
		else if (utile_per_contare == "S") N_tipi[_SOLENOID_]++;
		else if (utile_per_contare == "F") N_tipi[_FOCUSING_]++;
		else if (utile_per_contare == "D") N_tipi[_DEFOCUSING_]++;
		else log_file << "Something is wrong in the lattice file: " << utile_per_contare << " is not recognized\n";
    }


	Drift * drift = new Drift[N_tipi[_DRIFT_]];
	Solenoid * solenoid = new Solenoid[N_tipi[_SOLENOID_]];
	Focusing * focusing = new Focusing[N_tipi[_FOCUSING_]];
	Defocusing * defocusing = new Defocusing[N_tipi[_DEFOCUSING_]];

	Elemento_magnetico ** elementi = new Elemento_magnetico * [N_TIPI_MAGNETICI];  // puntatori alla classe base
	elementi[_DRIFT_] = drift;
	elementi[_SOLENOID_] = solenoid;
	elementi[_FOCUSING_] = focusing;
	elementi[_DEFOCUSING_] = defocusing;

	parametri_lattice.clear(), parametri_lattice.seekg(0, std::ios::beg);

	int numero_elementi_lattice= 0;
	for(int i=0; i < N_TIPI_MAGNETICI; i++) numero_elementi_lattice += N_tipi[i];

	if (numero_elementi_lattice == 0)
	{
		std::cerr << "Lattice must not be empty!" << std::endl;
		return -7;
	}

	double *param = new double[N_PARAMETRI_LATTICINO*numero_elementi_lattice];
	double lunghezza_totale = 0.0;
	double temp;
	int contatore = 0;
	while(contatore < numero_elementi_lattice)
	{
		parametri_lattice >> utile_per_contare;
		parametri_lattice >> parametro_1;
		parametri_lattice >> parametro_2;
		parametro_1s = parametro_1.c_str();
		parametro_2s = parametro_2.c_str();
		parametro_letto_1 = atof(parametro_1s);
		parametro_letto_2 = atof(parametro_2s);


		if (utile_per_contare == "O")
		{
			drift -> setvalues(parametro_letto_1, parametro_letto_2, lunghezza_totale);
			param[N_PARAMETRI_LATTICINO*contatore] = (double) _DRIFT_;		// tipo
			param[N_PARAMETRI_LATTICINO*contatore+1] = lunghezza_totale;	// posizione iniziale
			param[N_PARAMETRI_LATTICINO*contatore+2] = parametro_letto_1;	// lunghezza
			param[N_PARAMETRI_LATTICINO*contatore+3] = 0.0;					// campo
			drift++;
			lunghezza_totale += parametro_letto_1;
		}
		else if (utile_per_contare == "S")
		{
			solenoid -> setvalues(parametro_letto_1, parametro_letto_2, lunghezza_totale);
			param[N_PARAMETRI_LATTICINO*contatore] = (double) _SOLENOID_;	// tipo
			param[N_PARAMETRI_LATTICINO*contatore+1] = lunghezza_totale;	// posizione iniziale
			param[N_PARAMETRI_LATTICINO*contatore+2] = parametro_letto_1;	// lunghezza
			temp = CHARGE * parametro_letto_2 * FROM_TESLA_TO_GAUSS / (MP*C*C);
			param[N_PARAMETRI_LATTICINO*contatore+3] = temp;				// campo
			solenoid++;
			lunghezza_totale += parametro_letto_1;
		}
		else if (utile_per_contare == "F")
		{
			focusing -> setvalues(parametro_letto_1, parametro_letto_2, lunghezza_totale);
			param[N_PARAMETRI_LATTICINO*contatore] = (double) _FOCUSING_;	// tipo
			param[N_PARAMETRI_LATTICINO*contatore+1] = lunghezza_totale;	// posizione iniziale
			param[N_PARAMETRI_LATTICINO*contatore+2] = parametro_letto_1;	// lunghezza
			temp = parametro_letto_2 * FROM_TESLA_TO_GAUSS / FROM_M_TO_CM;
			if (parametro_letto_1 > sqrt(temp)) log_file << "Attention: you're outside of the range of the thin lens approximation" << std::endl;
			temp = CHARGE * parametro_letto_2 * FROM_TESLA_TO_GAUSS / (FROM_M_TO_CM * MP * C * C);
			param[N_PARAMETRI_LATTICINO*contatore+3] = temp;				// campo
			focusing++;
			lunghezza_totale += parametro_letto_1;
		}
		else if (utile_per_contare == "D")
		{
			defocusing -> setvalues(parametro_letto_1, parametro_letto_2, lunghezza_totale);
			param[N_PARAMETRI_LATTICINO*contatore] = (double) _DEFOCUSING_;	// tipo
			param[N_PARAMETRI_LATTICINO*contatore+1] = lunghezza_totale;	// posizione iniziale
			param[N_PARAMETRI_LATTICINO*contatore+2] = parametro_letto_1;	// lunghezza
			temp = parametro_letto_2 * FROM_TESLA_TO_GAUSS / FROM_M_TO_CM;
			if (parametro_letto_1 > sqrt(temp)) log_file << "Attention: you're outside of the range of the thin lens approximation" << std::endl;
			temp = CHARGE * parametro_letto_2 * FROM_TESLA_TO_GAUSS / (FROM_M_TO_CM * MP * C * C);
			param[N_PARAMETRI_LATTICINO*contatore+3] = temp;				// campo
			defocusing++;
			lunghezza_totale += parametro_letto_1;
		}
		else log_file << "Something is wrong in the lattice file: " << utile_per_contare << " is not recognized\n";
		contatore++;
	}

	create_gnuplot_file(gnuplot_filename_string, run_name, param);

	double beta = sqrt(2.0*E0/MP_MEV);
	double pz0 = sqrt(beta*beta - px0*px0 - py0*py0);

	double * x = new double[N_DIMENSIONI_SPAZIO_FASI+1];			// x, x', y, y', z, pz, beta_tot
	x[0] = x0;
	x[1] = px0/pz0;
	x[2] = y0;
	x[3] = py0/pz0;
	x[4] = z0;
	x[5] = pz0;
	x[6] = beta;
	output_file << x[0] << "\t" << x[2] << "\t" << x[4] << "\t" << x[1]*x[5] << "\t" << x[3]*x[5] << "\t" << x[5] << std::endl;

	double * xstep = new double[(N_DIMENSIONI_SPAZIO_FASI+1)*numero_elementi_lattice];

	double * xnew = new double[(N_DIMENSIONI_SPAZIO_FASI+1)];
	for (int i = 0; i < N_DIMENSIONI_SPAZIO_FASI; i++) xnew[i] = 0.0;

	int numero_step_per_elemento = nsteps / numero_elementi_lattice;
	numero_step_per_elemento++;

	double posizione = 0.0;
	double deltaZ;

	for (int i = 0; i < numero_elementi_lattice; i++)
	{
		deltaZ = param[N_PARAMETRI_LATTICINO*i+2] / (double) numero_step_per_elemento;
		for (int j = 0; j < numero_step_per_elemento; j++)
		{
			posizione += deltaZ;
			if (param[N_PARAMETRI_LATTICINO*i] < _DRIFT_CHECK_) Drift::mappa_(x, xnew, posizione,  param+N_PARAMETRI_LATTICINO*i);
			else if (param[N_PARAMETRI_LATTICINO*i] < _SOLENOID_CHECK_) Solenoid::mappa_(x, xnew, posizione,  param+N_PARAMETRI_LATTICINO*i);
			else if (param[N_PARAMETRI_LATTICINO*i] < _FOCUSING_CHECK_) Focusing::mappa_(x, xnew, posizione,  param+N_PARAMETRI_LATTICINO*i);
			else if (param[N_PARAMETRI_LATTICINO*i] < _DEFOCUSING_CHECK_) Defocusing::mappa_(x, xnew, posizione,  param+N_PARAMETRI_LATTICINO*i);
			else std::cerr << "Out of lattice" << std::endl;
			output_file << x[0] << "\t" << x[2] << "\t" << x[4] << "\t" << x[1]*x[5] << "\t" << x[3]*x[5] << "\t" << x[5] << std::endl;
		}
		for (int j = 0; j < N_DIMENSIONI_SPAZIO_FASI; j++) xstep[j+i] = xnew[j];
		xstep[6+i] = x[6];
		output_file << xnew[0] << "\t" << xnew[2] << "\t" << xnew[4] << "\t" << xnew[1]*xnew[5] << "\t" << xnew[3]*xnew[5] << "\t" << xnew[5] << std::endl;
	}

	log_file.close();
	output_file.close();
}

