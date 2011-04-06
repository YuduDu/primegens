/*Calculate energy for cross-hybridizations.*/

/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

    This file is part of the oligotm library.

    The oligotm library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    The oligotm library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the oligtm library (file gpl.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <limits.h>
#include <math.h>
#include <string.h>
#include "oligotm.h"

#define A_CHAR 'A'
#define G_CHAR 'G'
#define T_CHAR 'T'
#define C_CHAR 'C'
#define N_CHAR 'N'

#define H_CHAR 'H'
#define I_CHAR 'I'
#define J_CHAR 'J'
#define K_CHAR 'K'
#define L_CHAR 'L'
#define M_CHAR 'M'
#define V_CHAR 'V'
#define O_CHAR 'O'
#define P_CHAR 'P'
#define Q_CHAR 'Q'
#define R_CHAR 'R'
#define S_CHAR 'S'


#define CATID5(A,B,C,D,E) A##B##C##D##E
#define CATID2(A,B) A##B

#define STATE(LAST)     \
   CATID2(LAST,_STATE): \
   c = *s; s++;         \
   DO_PAIR(LAST,A)      \
   else DO_PAIR(LAST,T) \
   else DO_PAIR(LAST,G) \
   else DO_PAIR(LAST,C) \
   else DO_PAIR(LAST,N) \
   else if ('\0' == c)  \
             goto DONE; \
   else goto ERROR \

#define STATE2(LAST)     \
   CATID2(LAST,_STATE2): \
   c = *s; s++;         \
   DO_PAIR2(LAST,A)      \
   else DO_PAIR2(LAST,T) \
   else DO_PAIR2(LAST,G) \
   else DO_PAIR2(LAST,C) \
   else DO_PAIR2(LAST,N) \
   else DO_PAIR2(LAST,H) \
   else DO_PAIR2(LAST,I) \
   else DO_PAIR2(LAST,J) \
   else DO_PAIR2(LAST,K) \
   else DO_PAIR2(LAST,L) \
   else DO_PAIR2(LAST,M) \
   else DO_PAIR2(LAST,N) \
   else DO_PAIR2(LAST,O) \
   else DO_PAIR2(LAST,P) \
   else DO_PAIR2(LAST,Q) \
   else DO_PAIR2(LAST,R) \
   else DO_PAIR2(LAST,S) \
   else if ('\0' == c)  \
   goto DONE; \
   else goto ERROR \

/*
 * Two tables of nearest-neighbor parameters for di-nucleotide
 * base pairs.
 *
 * These are included in this file because they are not needed by
 * clients (callers) of oligtm().
 */

/* Table 1 (old parameters):
 * See table 2 in the paper [Breslauer KJ, Frank R, Blöcker H and
 * Marky LA (1986) "Predicting DNA duplex stability from the base
 * sequence" Proc Natl Acad Sci 83:4746-50
 * http://dx.doi.org/10.1073/pnas.83.11.3746]
 */

/* Delta G's of disruption * 1000. */
#define G_A_A  1900
#define G_A_C  1300
#define G_A_G  1600
#define G_A_T  1500
#define G_A_N  1575

#define G_C_A  1900
#define G_C_C  3100
#define G_C_G  3600
#define G_C_T  1600
#define G_C_N  2550

#define G_G_A  1600
#define G_G_C  3100
#define G_G_G  3100
#define G_G_T  1300
#define G_G_N  2275

#define G_T_A   900
#define G_T_C  1600
#define G_T_G  1900
#define G_T_T  1900
#define G_T_N  1575

#define G_N_A  1575
#define G_N_C  2275
#define G_N_G  2550
#define G_N_T  1575
#define G_N_N  1994

/* Table 2, new parameters:
 * Tables of nearest-neighbor thermodynamics for DNA bases, from the
 * paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
 * and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
 * Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]
 */

/* Delta G's of disruption * 1000. */
#define DG_A_A  1000
#define DG_A_C  1440
#define DG_A_G  1280
#define DG_A_T  880
#define DG_A_N  880

#define DG_C_A  1450
#define DG_C_C  1840
#define DG_C_G  2170
#define DG_C_T  1280
#define DG_C_N  1450

#define DG_G_A  1300
#define DG_G_C  2240
#define DG_G_G  1840
#define DG_G_T  1440
#define DG_G_N  1300

#define DG_T_A   580
#define DG_T_C  1300
#define DG_T_G  1450
#define DG_T_T  1000
#define DG_T_N   580

#define DG_N_A   580
#define DG_N_C  1300
#define DG_N_G  1280
#define DG_N_T   880
#define DG_N_N   580

/* End of tables nearest-neighbor parameter. */


#define DG_H_T -610
#define DG_H_G -430
#define DG_H_C -170
#define DG_H_A -690
#define DG_H_N -475
#define DG_A_H -610
#define DG_C_H -430
#define DG_G_H -170
#define DG_T_H -690
#define DG_N_H -475

#define DG_A_I -880
#define DG_C_I -750
#define DG_G_I -810
#define DG_T_I -920
#define DG_N_I -840
#define DG_I_T -770
#define DG_I_G -790
#define DG_I_C -470
#define DG_I_A -1330
#define DG_I_N -940

#define DG_A_J -140
#define DG_C_J -30
#define DG_G_J 250
#define DG_T_J -420
#define DG_N_J -85
#define DG_J_T -20
#define DG_J_G -110
#define DG_J_C 520
#define DG_J_A -740
#define DG_J_N -87.5

#define DG_K_T -1330
#define DG_K_G -700
#define DG_K_C -790
#define DG_K_A -1050
#define DG_K_N -967.5
#define DG_A_K -1330
#define DG_C_K -700
#define DG_G_K -790
#define DG_T_K -1050
#define DG_N_K -967.5

#define DG_L_T -880
#define DG_L_G -750
#define DG_L_C -810
#define DG_L_A -920
#define DG_L_N -840
#define DG_A_L -770
#define DG_C_L -790
#define DG_G_L -470
#define DG_T_L -1330
#define DG_N_L -940

#define DG_M_T -730
#define DG_M_G -400
#define DG_M_C -980
#define DG_M_A -750
#define DG_M_N -715
#define DG_A_M -640
#define DG_C_M -620
#define DG_G_M -620
#define DG_T_M -970
#define DG_N_M -712.5

#define DG_V_T 130
#define DG_V_G 110
#define DG_V_C 1110
#define DG_V_A -440
#define DG_V_N 227.5
#define DG_A_V 130
#define DG_C_V 110
#define DG_G_V 1110
#define DG_T_V -440
#define DG_N_V 227.5

#define DG_A_O -20
#define DG_C_O -110
#define DG_G_O 520
#define DG_T_O -740
#define DG_N_O -87.5
#define DG_O_T -140
#define DG_O_G -30
#define DG_O_C 250
#define DG_O_A -420
#define DG_O_N -85

#define DG_P_T -70
#define DG_P_G 320
#define DG_P_C 590
#define DG_P_A -340
#define DG_P_N 125
#define DG_A_P -710
#define DG_C_P 470
#define DG_G_P -80
#define DG_T_P -430
#define DG_N_P -187.5

#define DG_Q_T -690
#define DG_Q_G 120
#define DG_Q_C -450
#define DG_Q_A -680
#define DG_Q_N -425
#define DG_A_Q -690
#define DG_C_Q 0120
#define DG_G_Q -450
#define DG_T_Q -680
#define DG_N_Q -425

#define DG_R_T -640
#define DG_R_G -620
#define DG_R_C -620
#define DG_R_A -970
#define DG_R_N -712.5
#define DG_A_R -730
#define DG_C_R -400
#define DG_G_R -980
#define DG_T_R -750
#define DG_N_R -715

#define DG_S_T -710
#define DG_S_G 470
#define DG_S_C -80
#define DG_S_A -430
#define DG_S_N -187.5
#define DG_A_S -70
#define DG_C_S 320
#define DG_G_S 590
#define DG_T_S -340
#define DG_N_S 125

#define DG_H_S 0
#define DG_I_S 0
#define DG_J_S 0
#define DG_K_S 0
#define DG_L_S 0
#define DG_M_S 0
#define DG_V_S 0
#define DG_O_S 0
#define DG_P_S 0
#define DG_Q_S 0
#define DG_R_S 0
#define DG_S_S 0


#define DG_H_H 0
#define DG_H_I 0
#define DG_H_J 0
#define DG_H_K 0
#define DG_H_L 0
#define DG_H_M 0
#define DG_H_V 0
#define DG_H_O 0
#define DG_H_P 0
#define DG_H_Q 0
#define DG_H_R 0
#define DG_H_S 0

#define DG_I_H 0
#define DG_I_I 0
#define DG_I_J 0
#define DG_I_K 0
#define DG_I_L 0
#define DG_I_M 0
#define DG_I_V 0
#define DG_I_O 0
#define DG_I_P 0
#define DG_I_Q 0
#define DG_I_R 0
#define DG_I_S 0

#define DG_J_H 0
#define DG_J_I 0
#define DG_J_J 0
#define DG_J_K 0
#define DG_J_L 0
#define DG_J_M 0
#define DG_J_V 0
#define DG_J_O 0
#define DG_J_P 0
#define DG_J_Q 0
#define DG_J_R 0
#define DG_J_S 0

#define DG_K_H 0
#define DG_K_I 0
#define DG_K_J 0
#define DG_K_K 0
#define DG_K_L 0
#define DG_K_M 0
#define DG_K_V 0
#define DG_K_O 0
#define DG_K_P 0
#define DG_K_Q 0
#define DG_K_R 0
#define DG_K_S 0

#define DG_L_H 0
#define DG_L_I 0
#define DG_L_J 0
#define DG_L_K 0
#define DG_L_L 0
#define DG_L_M 0
#define DG_L_V 0
#define DG_L_O 0
#define DG_L_P 0
#define DG_L_Q 0
#define DG_L_R 0
#define DG_L_S 0
 
#define DG_M_H 0
#define DG_M_I 0
#define DG_M_J 0
#define DG_M_K 0
#define DG_M_L 0
#define DG_M_M 0
#define DG_M_V 0
#define DG_M_O 0
#define DG_M_P 0
#define DG_M_Q 0
#define DG_M_R 0
#define DG_M_S 0
 
#define DG_V_H 0
#define DG_V_I 0
#define DG_V_J 0
#define DG_V_K 0
#define DG_V_L 0
#define DG_V_M 0
#define DG_V_V 0
#define DG_V_O 0
#define DG_V_P 0
#define DG_V_Q 0
#define DG_V_R 0
#define DG_V_S 0

#define DG_O_H 0
#define DG_O_I 0
#define DG_O_J 0
#define DG_O_K 0
#define DG_O_L 0
#define DG_O_M 0
#define DG_O_V 0
#define DG_O_O 0
#define DG_O_P 0
#define DG_O_Q 0
#define DG_O_R 0
#define DG_O_S 0

#define DG_P_H 0
#define DG_P_I 0
#define DG_P_J 0
#define DG_P_K 0
#define DG_P_L 0
#define DG_P_M 0
#define DG_P_V 0
#define DG_P_O 0
#define DG_P_P 0
#define DG_P_Q 0
#define DG_P_R 0
#define DG_P_S 0

#define DG_Q_H 0
#define DG_Q_I 0
#define DG_Q_J 0
#define DG_Q_K 0
#define DG_Q_L 0
#define DG_Q_M 0
#define DG_Q_V 0
#define DG_Q_O 0
#define DG_Q_P 0
#define DG_Q_Q 0
#define DG_Q_R 0
#define DG_Q_S 0

#define DG_R_H 0
#define DG_R_I 0
#define DG_R_J 0
#define DG_R_K 0
#define DG_R_L 0
#define DG_R_M 0
#define DG_R_V 0
#define DG_R_O 0
#define DG_R_P 0
#define DG_R_Q 0
#define DG_R_R 0
#define DG_R_S 0

#define DG_S_H 0
#define DG_S_I 0
#define DG_S_J 0
#define DG_S_K 0
#define DG_S_L 0
#define DG_S_M 0
#define DG_S_V 0
#define DG_S_O 0
#define DG_S_P 0
#define DG_S_Q 0
#define DG_S_R 0
#define DG_S_S 0


/* Calculate the melting temperature of oligo s.  See
   oligotm.h for documentation of arguments.
*/

#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dg += CATID5(G,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

/*
#define DO_PAIR2(LAST,THIS)          \
     if (CATID2(THIS,_CHAR) == c) { \
     dg += CATID5(DG,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE2);      \
}
*/

/*
#define DO_PAIR2(LAST,THIS)          \
     if (CATID2(THIS,_CHAR) == c) { \
     if (LAST=='S'){ dg += CATID5(DG,_,LAST_OF_LAST,_,P);}\
     else { dg += CATID5(DG,_,LAST,_,THIS);}\
     LAST_OF_LAST = LAST;            \
    goto CATID2(THIS,_STATE2);      \
}
*/

#define DO_PAIR3(LAST_OF_LAST,OTHER) \
       if(CATID2(A,_CHAR)==LAST_OF_LAST) dg += CATID5(DG,_,A,_,OTHER); \
       if(CATID2(T,_CHAR)==LAST_OF_LAST) dg += CATID5(DG,_,T,_,OTHER); \
       if(CATID2(G,_CHAR)==LAST_OF_LAST) dg += CATID5(DG,_,G,_,OTHER); \
       if(CATID2(C,_CHAR)==LAST_OF_LAST) dg += CATID5(DG,_,C,_,OTHER); \
       if(CATID2(N,_CHAR)==LAST_OF_LAST) dg += CATID5(DG,_,N,_,OTHER); \
}

#define DO_PAIR2(LAST,THIS)          \
       if(CATID2(THIS,_CHAR) == c) { \
	 if(CATID2(LAST,_CHAR) == CATID2(H,_CHAR)) {\
	    DO_PAIR3(LAST_OF_LAST,H)\
         if(CATID2(LAST,_CHAR) == CATID2(I,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,L)\
         if(CATID2(LAST,_CHAR) == CATID2(J,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,O)\
         if(CATID2(LAST,_CHAR) == CATID2(K,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,K)\
         if(CATID2(LAST,_CHAR) == CATID2(L,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,I)\
         if(CATID2(LAST,_CHAR) == CATID2(M,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,R)\
         if(CATID2(LAST,_CHAR) == CATID2(V,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,V)\
         if(CATID2(LAST,_CHAR) == CATID2(O,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,J)\
         if(CATID2(LAST,_CHAR) == CATID2(P,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,S)\
         if(CATID2(LAST,_CHAR) == CATID2(Q,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,Q)\
         if(CATID2(LAST,_CHAR) == CATID2(R,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,M)\
         if(CATID2(LAST,_CHAR) == CATID2(S,_CHAR)) {\
            DO_PAIR3(LAST_OF_LAST,P)\
	 else dg += CATID5(DG,_,LAST,_,THIS); \
	 LAST_OF_LAST = CATID2(LAST,_CHAR);\
	 goto CATID2(THIS,_STATE2);      \
}


float 
oligodg(s, tm_santalucia)
    const char *s;       /* The sequence. */
    int tm_santalucia; 
{
   register float dg = 0;
   register char c;
   register char LAST_OF_LAST = 'A';

//printf("s=%s\n",s);
   /* Use a finite-state machine (DFA) to calucluate dg s. */
   c = *s; s++;
   if(tm_santalucia ==1) {      
      dg=-1960; /* Initial dG */
      if(c == 'A' || c == 'T')  {
	 dg += -50; /* terminal AT penalty */
      }
      if (c == 'A') goto A_STATE2;
      else if (c == 'G') goto G_STATE2;
      else if (c == 'T') goto T_STATE2;
      else if (c == 'C') goto C_STATE2;
      else if (c == 'N') goto N_STATE2;
      else if (c == 'H') goto S_STATE2;
      else if (c == 'I') goto S_STATE2;
      else if (c == 'J') goto S_STATE2;
      else if (c == 'K') goto S_STATE2;
      else if (c == 'L') goto S_STATE2;
      else if (c == 'M') goto S_STATE2;
      else if (c == 'V') goto S_STATE2;
      else if (c == 'O') goto S_STATE2;
      else if (c == 'P') goto S_STATE2;
      else if (c == 'Q') goto S_STATE2;
      else if (c == 'R') goto S_STATE2;
      else if (c == 'S') goto S_STATE2;
      else goto ERROR;
      STATE2(A);
      STATE2(T);
      STATE2(G);
      STATE2(C);
      STATE2(N);
      STATE2(H);
      STATE2(I);
      STATE2(J);
      STATE2(K);
      STATE2(L);
      STATE2(M);
      STATE2(V);
      STATE2(O);
      STATE2(P);
      STATE2(Q);
      STATE2(R);
      STATE2(S);
     } else {
    if (c == 'A') goto A_STATE;
    else if (c == 'G') goto G_STATE;
    else if (c == 'T') goto T_STATE;
    else if (c == 'C') goto C_STATE;
    else if (c == 'N') goto N_STATE;
    else goto ERROR;
    STATE(A);
    STATE(T);
    STATE(G);
    STATE(C);
    STATE(N);

     }
DONE:  /* dg is now computed for the given sequence. */
   if(tm_santalucia ==1) {
      int sym;
      --s; --s; c = *s;
      if(c == 'A' || c == 'T')  {
	 dg += -50; /* terminal AT penalty */
      }
      sym = symmetry(s);
      if(sym==1)   {
	 dg +=-430; /* symmetry correction for dG */
      }
   }
//printf("dg=%f\n",dg);
   return dg/1000.0;

 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}

double end_oligodg(s, len, tm_santalucia)
  const char *s;  
  int len; /* The number of characters to return. */
  int tm_santalucia;
{
  int x = strlen(s);

  if (tm_santalucia != TM_METHOD_BRESLAUER
      && tm_santalucia != TM_METHOD_SANTALUCIA)
    return OLIGOTM_ERROR;

  return 
    x < len 
    ? oligodg(s,tm_santalucia) :
    oligodg(s + (x - len),tm_santalucia);
}

 /* Return 1 if string is symmetrical, 0 otherwise. */ 
int symmetry(const char* seq) {
   register char s;
   register char e;
   const char *seq_end=seq;
   int i = 0;
   int seq_len=strlen(seq);
   int mp = seq_len/2;
   if(seq_len%2==1) {
      return 0;
   }
   seq_end+=seq_len;
   seq_end--;
   while(i<mp) {
      i++;
      s=*seq;
      e=*seq_end;
      if ((s=='A' && e!='T') 
	  || (s=='T' && e!='A') 
	  || (e=='A' && s!='T') 
	  || (e=='T' && s!='A')) {
	 return 0;
      }
      if ((s=='C' && e!='G')
	  || (s=='G' && e!='C')
	  || (e=='C' && s!='G')
	  || (e=='G' && s!='C')) {
	 return 0;
      }
      seq++;
      seq_end--;
   }
   return 1;
}


