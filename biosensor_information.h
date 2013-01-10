#ifndef BIOSENSORINFORMATION_H
#define BIOSENSORINFORMATION_H

enum resp_method
{
	DEFAULT_TIME, //0 - iki pusiausvyros
	MIN_TIME,     //1 - iki pusiausvyros su nurodytu minimaliu laiku
	FIXED_TIME    //2 - fiksuotas laikas
};

struct layer_params 
{
	//Laukas nurodo ar tai fermento sluoksnis
	int enz_layer;
	//Difuzijos koeficientai (cm^2/s)
	double Ds;
	double Dp;
	//Sluoksnio storis (cm)
	double d;
	//Fermento koncentracija (mol/cm^3)
	double e0;
};

struct bio_params 
{
	//Laukas nurodo ar naudojama išreikštinė schema (priešingu atveju bus naudojama neišreikštinė schema)
	int explicit_scheme;
	//Laukas nurodo ar vyksta substrato inhibicija
	int subs_inh;
	//Laukas nurodo ar vyksta produkto inhibicija
	int prod_inh;
	//Reakcijos greičio konstanta k2 (s^-1)
	double k2;
	//Pusiausvyros konstantos (mol/cm^3)
	double km;
	double ks;
	double kp;
	//Žingsnis pagal laiką (s)
	double dt;
	//Į kiek dalių dalinami sluoksniai
	int n;
	//Metodas, kuriuo bus nustatomas atsako laikas:
	enum resp_method resp_t_meth;
	//Minimalus atsako laikas (s)
	double min_t;
	//Fiksuotas atsako laikas (s)
	double resp_t;
	//Išvedimo failas
	char *out_file_name;
	//Elektronų, dalyvaujančių krūvio pernešime, skaičius
	int ne;
	//Substrato koncentracija tirpale (mol/cm^3)
	double s0;
	//Produkto koncentracija tirpale (mol/cm^3)
	double p0;
	//Biojutiklio sluoksnių skaičius
	int layer_count;
	//Biojutiklio sluoksnių masyvas
	struct layer_params *layers;
};

#endif // BIOSENSORINFORMATION_H
