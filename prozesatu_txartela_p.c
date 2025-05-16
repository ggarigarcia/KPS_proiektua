/*
    prozesatu_txartela_s.c    (txip_banaketa)  

    txartela prozesatzeko funtzio nagusiak

******************************************************************/

#include <mpi.h>

#include "definizioak.h"


extern int pid, prk, err_kop;
extern int *sendcounts, *displs;


/*****************************************/
/**
 * @brief Beroa txertatu txipen posizioetan azpi_sareta matrizean, azpi_sareta matrizearen erdialdea freskatu
 */
void txertatu_beroa (struct info_param param, float *sareta, float *sareta_txip)
{
  int  i, j, a, b;


  // txertatu beroa txipen posizioetan, tenperatura-diferentziaren arabera
  for (i= 0; i< err_kop; i++)
  for (j=1; j<ZUT-1; j++) 
    if (sareta_txip[i*ZUT+j] > sareta[(i+1)*ZUT+j])
      sareta[(i+1)*ZUT+j] += 0.05 * (sareta_txip[i*ZUT+j] - sareta[(i+1)*ZUT+j]);


  // freskatu txartelaren erdialdea (haizegailuak erdiko zerrenda bertikalean) 
  a = 0.44*(ZUT-2) + 1;
  b = 0.56*(ZUT-2) + 1;

  for (i=1; i <= err_kop; i++)
    for (j=a; j<b; j++)
      sareta[i*ZUT+j] -= 0.01 * (sareta[i*ZUT+j] - param.t_kanpo);
}


/**
 * @brief Prozesuen artean halo/ghost errenkadak elkartrukatu, funtzio SINKRONOAK erabiliz
 */
void errenkadak_elkartrukatu_ssend(struct info_param param, float *sareta) {

  // * Bikoitiek SSEND (azken errenkada PROPIOA), bakoitiek RECV (lehen errenkada)
    // * Bakoitiek SSEND (lehenengo errenkada PROPIOA), bikoitiek RECV (azken errenkada)

  // * Bikoitiek SSEND (lehen errenkada PROPIOA), bakoitiek RECV (azken errenkada)
    // * Bakoitiek SSEND (azken errenkada PROPIOA), bikoitiek RECV (lehen errenkada)

  if(pid % 2 == 0) {

    if(pid != prk-1) {
      MPI_Ssend(&sareta[err_kop * ZUT], ZUT, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD);
      MPI_Recv(&sareta[(err_kop+1) * ZUT], ZUT, MPI_FLOAT, pid+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(pid != 0) {
      MPI_Ssend(&sareta[ZUT], ZUT, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD);
      MPI_Recv(&sareta[0], ZUT, MPI_FLOAT, pid-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  } else {

    if(pid != 0) {
      MPI_Recv(&sareta[0], ZUT, MPI_FLOAT, pid-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Ssend(&sareta[ZUT], ZUT, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD);
    }

    if(pid != prk-1) {
      MPI_Recv(&sareta[(err_kop+1) * ZUT], ZUT, MPI_FLOAT, pid+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Ssend(&sareta[err_kop * ZUT], ZUT, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD);
    }

  }
}


/**
 * @brief Prozesuen artean errenkadak elkartrukatu, funtzio ASINKRONOAK erabiliz
 */
void errenkadak_elkartrukatu_isend(struct info_param param, float *sareta) {
 
  // * Errenkaden zain jarri
  // * Errenkadak bidali
 
 
  MPI_Request reqs[4];
  int req_count = 0;

  
  if (pid != 0) {
      MPI_Irecv(&sareta[0], ZUT, MPI_FLOAT, pid - 1, 0, MPI_COMM_WORLD, &reqs[req_count++]); 
  }
  if (pid != prk - 1) {
      MPI_Irecv(&sareta[(err_kop + 1) * ZUT], ZUT, MPI_FLOAT, pid + 1, 1, MPI_COMM_WORLD, &reqs[req_count++]); 
  }

  
  if (pid != 0) {
      MPI_Isend(&sareta[ZUT], ZUT, MPI_FLOAT, pid - 1, 1, MPI_COMM_WORLD, &reqs[req_count++]); 
  }
  if (pid != prk - 1) {
      MPI_Isend(&sareta[err_kop * ZUT], ZUT, MPI_FLOAT, pid + 1, 0, MPI_COMM_WORLD, &reqs[req_count++]); 
  }

  
  MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
}



/*****************************************/
/**
 * @brief azpi_sareta matrizean bero-difusioaren simulazioa
 */
double barreiatu_beroa (struct info_param param, float *sareta, float *sareta_lag)
{
  int     i, j;
  float   T;
  double  Tosoa;


  // beroaren difusioa: puntu bakoitzaren tenperatura aldatzen da auzokideen tenperaturen arabera
  // aldaketa sareta_lag matrizean egiten da; bukaeran, datuak kopiatzen dira sareta matrizera

  Tosoa = 0.0;
  for (i=1; i <= err_kop; i++)
  for (j=1; j<ZUT-1; j++)
  {
    T = sareta[i*ZUT + j] + 0.10 * ( 
          sareta[(i-1)*ZUT + j-1] + sareta[(i-1)*ZUT + j] + sareta[(i-1)*ZUT + j+1] +
          sareta[    i*ZUT + j-1] + sareta[    i*ZUT + j+1] +
          sareta[(i+1)*ZUT + j-1] + sareta[(i+1)*ZUT + j] + sareta[(i+1)*ZUT + j+1] 
          - 8*sareta[i*ZUT + j] 
          );

    sareta_lag[i*ZUT + j] = T;
    Tosoa += T;
  }

  // saretarako balio berriak
  for (i=1; i <= err_kop; i++)
  for (j=1; j<ZUT-1; j++)
    sareta[i*ZUT+j] = sareta_lag[i*ZUT+j]; 
 
  return (Tosoa);
}



/*****************************************/
double kalkulatu_Tbb (struct info_param param, float *sareta, float *sareta_txip, float *sareta_lag, int elkartrukatu_mota)
{
  int     i, j, amaia, iterkop;
  double  Tosoa, Tbb, Tbb0 = param.t_kanpo;

  double Tosoa_global;

  // begizta nagusia, emaitzak konbergitu arte: t_delta edo iter_amax

  amaia = 0; 
  iterkop = 0;

  while (amaia == 0)
  {
    iterkop ++;

    // beroa eta haizeaz freskatu
    txertatu_beroa (param, sareta, sareta_txip);

    if (elkartrukatu_mota == 0) {
        errenkadak_elkartrukatu_ssend(param, sareta);
    } else {
        errenkadak_elkartrukatu_isend(param, sareta);
    }
 
    // beroa barreiatu eta batez besteko tenperatura lortu
    Tosoa = barreiatu_beroa (param, sareta, sareta_lag);

    // aztertu konbergentzia 10 iteraziotan behin
    if (iterkop % 10 == 0)
    {
      MPI_Allreduce(&Tosoa, &Tosoa_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      Tbb = Tosoa_global / ((ZUT-2)*(ERR-2));

      if ((fabs(Tbb - Tbb0) < param.t_delta) || (iterkop > param.iter_max))
           amaia = 1;
      else Tbb0 = Tbb;
    }
  } 

  if(pid == 0) printf ("  Iter_kop: %4d\t", iterkop);
  return (Tbb);
}

