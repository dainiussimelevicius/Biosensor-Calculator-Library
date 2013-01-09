#ifndef BIOSENSORINFORMATION_H
#define BIOSENSORINFORMATION_H

enum RESP_METHOD            
{
  DEFAULT_TIME, //0 - iki pusiausvyros
  MIN_TIME,     //1 - iki pusiausvyros su nurodytu minimaliu laiku
  FIXED_TIME    //2 - fiksuotas laikas
};

struct LayerInformation 
{
    //Laukas nurodo ar tai fermento sluoksnis
    int enzymeLayer;
    //Difuzijos koeficientai (cm^2/s)
    double Ds;
    double Dp;
    //Sluoksnio storis (cm)
    double d;
    //Fermento koncentracija (mol/cm^3)
    double e0;
};

struct BiosensorInformation 
{
    //Laukas nurodo ar naudojama išreikštinė schema (priešingu atveju bus naudojama neišreikštinė schema)
    int explicitScheme;
    //Laukas nurodo ar vyksta substrato inhibicija
    int substrateInhibition;
    //Laukas nurodo ar vyksta produkto inhibicija
    int productInhibition;
    //Reakcijos greičio konstanta k2 (s^-1)
    double k2;
    //Pusiausvyros konstantos (mol/cm^3)
    double kM;
    double kS;
    double kP;
    //Žingsnis pagal laiką (s)
    double timeStep;
    //Į kiek dalių dalinami sluoksniai
    int N;
    //Metodas, kuriuo bus nustatomas atsako laikas:
    enum RESP_METHOD responseTimeMethod;
    //Minimalus atsako laikas (s)
    double minTime;
    //Fiksuotas atsako laikas (s)
    double responseTime;
    //Išvedimo failas
    char *outputFileName;
    //Elektronų, dalyvaujančių krūvio pernešime, skaičius
    int ne;
    //Substrato koncentracija tirpale (mol/cm^3)
    double s0;
    //Produkto koncentracija tirpale (mol/cm^3)
    double p0;
    //Biojutiklio sluoksnių skaičius
    int noOfBiosensorLayers;
    //Biojutiklio sluoksnių masyvas
    struct LayerInformation *biosensorLayers;
};

#endif // BIOSENSORINFORMATION_H
