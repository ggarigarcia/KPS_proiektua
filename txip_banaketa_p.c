/*  KPS, Konputazio Paraleloko Sistemak -- IIG (Konputagailuen Ingeniaritza)
    MPI laborategia

    txip_banaketa_s.c 

    Beroaren difusioa 2 dimentsiotan      Serieko bertsioa

    Txartel elektroniko baten txipen posizioak analizatzen dira, sistemaren tenperatura minimoa lortzeko.
    Poisson motako metodoa erabiltzen da (sinplifikatuta), txartela hainbat puntutan diskretizatuta.

    Sarrera > txartela: txartelaren definizioa eta simulatu behar diren konfigurazioak
    Irteera > konfiguraziorik onena eta batez besteko tenperatura
       	      txartela_s.txipak: sarrerako egoera termikoa
              txartela_s.emaitzak: bukaerako egoera termikoa 
 
    Definizioak.h: zenbait aldagai eta datu-egituren definizioak

    Konpilatu bi fitxategi hauekin: 
      prozesatu_txartela_s.c: beroa txertatu, barreiatu, eta Tbb kalkulatu tenperatura egonkortu arte
      flag_s.c: hainbat funtzio laguntzaile

************************************************************************************************/


#include <stdio.h>
#include <values.h>
#include <time.h>
#include <mpi.h>

#include "definizioak.h"
#include "flag_s.h"
#include "prozesatu_txartela_s.h"


int pid, prk, err_kop;
int *sendcounts, *displs;
float *azpi_sareta_txip, *azpi_sareta, *azpi_sareta_lag;


/**
 * @brief Paralelizaziorako beharrezkoak diren aldagaiak prestatu
 */
void paralelizazioa_prestatu(struct info_param param) {

  // * sendcounts, displs (scatterv/gatherv)
  sendcounts = malloc(prk * sizeof(int));
  displs     = malloc(prk * sizeof(int));

  int zatidura = (ERR-2) / prk;
  int hondarra = (ERR-2) % prk;
  int offset   = 0;

  int i;
  for(i = 0; i < prk; i++) {

    sendcounts[i] = zatidura * ZUT;
    if(i < hondarra) sendcounts[i] += ZUT;
    displs[i] = offset;
    offset += sendcounts[i]; 
  }

  // * azpi saretak
  err_kop = sendcounts[pid] / ZUT; // ! elkartruke errenkadak EZ dira kontuan hartzen

  azpi_sareta_txip = malloc(err_kop * ZUT * sizeof(float));
    azpi_sareta      = malloc((err_kop+2) * ZUT * sizeof(float));
    azpi_sareta_lag  = malloc((err_kop+2) * ZUT * sizeof(float));
}

/*****************************************/
void hasi_sareta_txip (int konf, struct info_param param, struct info_txip *txipak, int **txip_koord, 
                       float *sareta_txip)
{
  int i, j, n;


  // hasierako tenperatura
  for (i=0; i<ERR; i++)
  for (j=0; j<ZUT; j++)  
    sareta_txip[i*ZUT+j] = param.t_kanpo;

  // txip bakoitzaren tenperatura, okupatzen dituzten posizioetan
  for (n=0; n<param.txip_kop; n++)
  {
    for (i = txip_koord[konf][2*n]   * param.eskala; i < (txip_koord[konf][2*n]   + txipak[n].altu)  * param.eskala; i++)
    for (j = txip_koord[konf][2*n+1] * param.eskala; j < (txip_koord[konf][2*n+1] + txipak[n].zabal) * param.eskala; j++)
      sareta_txip[(i+1)*ZUT+(j+1)] = txipak[n].t_txip;
  }
}  



/*****************************************/
/**
 * @brief Azpi saretak hasieratu, uneko konfigurazioaren hasierako tenperaturarekin
 */
void hasi_saretak (struct info_param param)
{
  // sareta eta sareta_lag datu-egiturak hasieratzen dira: kanpo tenperatura

  int  i, j;

  // * azpi_saretak
  for (i=0; i< err_kop+2; i++)
    for (j=0; j<ZUT; j++) 
    {
      azpi_sareta[i*ZUT+j] = param.t_kanpo;
      azpi_sareta_lag[i*ZUT+j] = param.t_kanpo;
    }
}



/*****************************************/
/*****************************************/
int main (int argc, char *argv[])
{
  struct info_param  param;
  struct info_txip   *txipak;
  int  		     **txip_koord;

  float  *sareta_txip, *sareta, *sareta_lag;
  struct info_emaitzak  ONENA;

  int     konf, i;
  double  Tbb, tex;
  struct timespec  t0, t1;
  int elkartrukatu_mota;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &prk);

  if (argc != 3) {
    if(pid == 0) printf ("\n\nERROREA: txartelaren fitxategia eta/edo elkartrukatze mota falta da (0: Ssend, 1: Isend) \n\n");
    MPI_Finalize();
    exit (-1);
  } 

  elkartrukatu_mota = atoi(argv[2]);
  if (elkartrukatu_mota != 0 && elkartrukatu_mota != 1) {
    if(pid == 0) printf ("\n\nERROREA: elkartrukatze mota baliogabea. Erabili 0 (Ssend) edo 1 (Isend).\n\n");
    MPI_Finalize();
    exit(-1);
  }

  // irakurri sarrera-datuak

  if(pid == 0) irakurri_datuak (argv[1], &param, &txipak, &txip_koord);
  MPI_Bcast(&param, sizeof(struct info_param), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(pid == 0) {
    printf ("\n===================================================================");
    printf ("\n  Beroa barreiatzea - Bertsio PARALELOA (%d prozesu) ", prk);
    printf ("\n  %d x %d puntu, %d txip", EMIN*param.eskala, ZMIN*param.eskala, param.txip_kop);
    printf ("\n  T_kanpo: %1.1f, Tmax_txip: %1.1f, T_delta: %1.3f, Iter_max: %d", param.t_kanpo,
          param.tmax_txip, param.t_delta, param.iter_max);
    printf ("\n===================================================================\n\n"); 

    clock_gettime (CLOCK_REALTIME, &t0);
  }

  if(pid == 0) {
    sareta_txip   = malloc (ERR*ZUT * sizeof(float));
    sareta        = malloc (ERR*ZUT * sizeof(float));
    //sareta_lag    = malloc (ERR*ZUT * sizeof(float));
    ONENA.bsareta = malloc (ERR*ZUT * sizeof(float));
    ONENA.csareta = malloc (ERR*ZUT * sizeof(float));
  }

  ONENA.Tbb = MAXDOUBLE;

  paralelizazioa_prestatu(param);

  // txipen konfigurazioak prozesatzeko begizta
  // ==========================================

  for (konf=0; konf<param.konf_kop; konf++)
  {
    // hasierako balioak 
    if(pid == 0) hasi_sareta_txip (konf, param, txipak, txip_koord, sareta_txip);

    MPI_Scatterv(&sareta_txip[ZUT], sendcounts, displs, MPI_FLOAT, &azpi_sareta_txip[0], sendcounts[pid], MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    hasi_saretak (param);

    // funtzio nagusia: beroa txertatu/barreiatu tenp. konbergitu arte (t_delta edo iter_max)
    Tbb = kalkulatu_Tbb (param, azpi_sareta, azpi_sareta_txip, azpi_sareta_lag, elkartrukatu_mota);
    if(pid == 0) printf ("  Konfig: %2d    Tbb: %1.2f\n", konf + 1, Tbb);
   
    // emaitza hobea bada, gorde
    if (Tbb < ONENA.Tbb) {
      
      MPI_Gatherv(&azpi_sareta[ZUT], sendcounts[pid], MPI_FLOAT, &sareta[ZUT], sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

      if(pid == 0) {
        konf_onena (konf, Tbb, param, sareta, sareta_txip, &ONENA);
      } else {
        ONENA.Tbb = Tbb;
        ONENA.konf = konf;
      }
    }
  }

  if(pid == 0) {
    clock_gettime (CLOCK_REALTIME, &t1);
    tex = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/(double)1e9;

    printf ("\n\n >>> Konfigurazio onena: %2d    Tbb: %1.2f\n", ONENA.konf + 1, ONENA.Tbb); 
    printf ("   > Tex (paralelo, %d prozesu): %1.3f s \n\n", prk, tex);
  }

  // idatzi fitxategian konfigurazio onenaren emaitzak
  if(pid == 0) idatzi_onena (param, &ONENA, argv[1]);

  if(pid == 0) {

    free (sareta_txip); 
    free (sareta); 
    //free (sareta_lag);
    free (ONENA.bsareta); 
    free (ONENA.csareta);
    free (txipak);
    for (i=0; i<param.konf_kop; i++) free (txip_koord[i]);
    free (txip_koord);
  }

  // free(azpi_sareta); free(azpi_sareta_lag); free(azpi_sareta_txip); free(sendcounts); free(displs); // ! SEGFAULT / MEMORY arazoak

  MPI_Finalize();
  return (0);
}

