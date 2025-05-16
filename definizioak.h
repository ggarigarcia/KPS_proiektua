/*
    definizioak.h  (txip_banaketa) 

******************************************************************/


// txartelaren oinarrizko tamaina
// errenkada eta zutabe kopurua txartelaren eskalaren arabera (1 - 12) 

#define EMIN 200
#define ZMIN 100

#define ERR (EMIN*param.eskala + 2)  // saretaren errenkada kopurua
#define ZUT (ZMIN*param.eskala + 2)  // saretaren zutabe kopurua


// datu-egiturak eta simulazioaren parametro orokorrak

// kanpoko tenperatura, txipen tenperatura maximoa, tenperatura-diferentzia txikiena
// konfigurazio kopurua, txip kopurua, iterazio kopuru maximoa, txartelaren eskala-faktorea

struct info_param 
{
  int    konf_kop, txip_kop, iter_max, eskala;
  float  t_kanpo, tmax_txip, t_delta;    
};


// txipen definizioa: tamaina eta tenperatura
 
struct info_txip 
{
  int    altu, zabal;
  float  t_txip;
};


// emaitzak: batez besteko tenp., azken sareta, hasierako sareta (txipak), konfigurazio zenbakia

struct info_emaitzak 
{
  double  Tbb;
  float   *bsareta;	// azken tenperatura banaketa
  float   *csareta;	// hasierako egoera
  int     konf;
};

