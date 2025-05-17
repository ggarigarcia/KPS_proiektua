/*
    flag_s.c

    txip_banaketa_s.c programarako funtzio laguntzaileak

**********************************************************************/


#include <stdio.h>
#include <values.h>

#include "definizioak.h"


/*****************************************/
void irakurri_datuak (char *fsar, struct info_param *param, struct info_txip **txipak, int ***txip_koord)
{
  int    i, j, altu, zabal;
  float  t_txip;
  FILE   *fd;


  fd = fopen (fsar, "r");
  if (fd == NULL) {
    printf ("\n\nERROREA: sarrera-fitxategia \n\n");
    exit (-1);
  }

  // irakurri simulazioaren parametroak
  fscanf (fd, "%d %d %d %f %f %f %d", &param->eskala, &param->konf_kop, &param->txip_kop, &param->t_kanpo, 
               &param->tmax_txip, &param->t_delta, &param->iter_max);

  if (param->eskala > 12) {
    printf ("\n\nERROREA: eskalaren balio maximoa, 12 \n\n");
    exit (-1);
  }

  // txipen tamainak eta tenperaturak
  *txipak = malloc (param->txip_kop * sizeof(struct info_txip));

  for (i=0; i<param->txip_kop; i++)
  {
    fscanf (fd, "%d %d %f", &altu, &zabal, &t_txip);
    (*txipak)[i].altu = altu;
    (*txipak)[i].zabal = zabal;
    (*txipak)[i].t_txip = t_txip;
  }

  // txipen posizioak
  *txip_koord = malloc (param->konf_kop * sizeof(int*));
  for (i=0; i<param->konf_kop; i++) 
    (*txip_koord)[i] = malloc (2 * param->txip_kop * sizeof(int));

  for (i=0; i<param->konf_kop; i++)
  for (j=0; j<param->txip_kop; j++)
    fscanf (fd, "%d %d", &(*txip_koord)[i][2*j], &(*txip_koord)[i][2*j+1]);

  fclose (fd);
}



/*****************************************/
void konf_onena (int konf, double Tbb, struct info_param param, 
                 float *sareta, float *sareta_txip, struct info_emaitzak *ONENA)
{
  // bsareta: tenperaturen azken banaketa; csareta: txipen konfigurazioa

  int  i, j;


  ONENA->Tbb = Tbb;
  ONENA->konf = konf;
  for (i=1; i<ERR-1; i++)
  for (j=1; j<ZUT-1; j++) 
  {
    ONENA->bsareta[i*ZUT+j] = sareta[i*ZUT+j];
    ONENA->csareta[i*ZUT+j] = sareta_txip[i*ZUT+j];
  }
}



/*****************************************/
void idatzi_sareta (FILE *fd, float *sareta, struct info_param param)
{
  int  i, j;


  // j - i ordena gero irudia horizontalki ikusteko
  for (j=ZUT-2; j>0; j--)
  {
    for (i=1; i<ERR-1; i++) fprintf (fd, "%1.2f ", sareta[i*ZUT+j]);
    fprintf (fd, "\n");
  }
  fprintf (fd, "\n");
}


/*****************************************/
void idatzi_onena (struct info_param param, struct info_emaitzak *ONENA, char *fsar)
{
  FILE  *fd;
  char  izena[100];


  // azken sareta
  sprintf (izena, "%s_p.emaitza", fsar);
  fd = fopen (izena, "w");
  fprintf (fd, "Tmin_hasi %1.1f  Tmax_hasi %1.1f  \n", param.t_kanpo, param.tmax_txip);
  fprintf (fd, "%d\t  %d \n", ZUT-2, ERR-2);

  idatzi_sareta (fd, ONENA->bsareta, param);

  fprintf (fd, "\n\n >>> Konfigurazio onena: %d\t Tbb: %1.2f\n\n", ONENA->konf+1, ONENA->Tbb);
  fclose (fd);

  // hasierako txip-sareta
  sprintf (izena, "%s_p.txipak", fsar);
  fd = fopen (izena, "w");
  fprintf (fd, "Tmin_txip %1.1f  Tmax_txip %1.1f  \n", param.t_kanpo, param.tmax_txip);
  fprintf (fd, "%d\t  %d \n", ZUT-2, ERR-2);

  idatzi_sareta (fd, ONENA->csareta, param);

  fclose (fd);
}
