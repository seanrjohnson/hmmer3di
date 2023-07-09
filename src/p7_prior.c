/* Mixture Dirichlet priors for profile HMMs.
 * 
 * SRE, Sat Mar 24 09:12:44 2007 [Janelia]
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"


/* Function:  p7_prior_CreateAmino()
 * Incept:    SRE, Sat Mar 24 09:35:36 2007 [Janelia]
 *
 * Purpose:   Creates the default mixture Dirichlet prior for protein
 *            sequences.
 *
 *            The transition priors (match, insert, delete) are all
 *            single Dirichlets, originally trained by Graeme
 *            Mitchison in the mid-1990's. Notes have been lost, but
 *            we believe they were trained on an early version of
 *            Pfam. 
 *            
 *            The match emission prior is a nine-component mixture
 *            from Kimmen Sjolander, who trained it on the Blocks9
 *            database \citep{Sjolander96}.
 *            
 *            The insert emission prior is a single Dirichlet with
 *            high $|\alpha|$, such that insert emission probabilities
 *            are essentially fixed by the prior, regardless of
 *            observed count data. The slightly polar parameterization
 *            was obtained by training on Pfam 1.0.
 *
 * Returns:   a pointer to the new <P7_PRIOR> structure.
 */
P7_PRIOR *
p7_prior_CreateAmino(void)
{
  P7_PRIOR *pri = NULL;
  int k;
  int status;
				/* default match mixture coefficients: [Sjolander96] */
        /* calculated by esl-mixdchlet from 3Di profiles from AlphaFold UniProt for Pfam 35.0, Sean Johnson June 4, 2023 */
  static double defmq[9] = {  0.1206, 0.0061, 0.1222, 0.0892, 0.2931, 0.1027, 0.0173, 0.1001, 0.1487 };

				/* default match mixture Dirichlet components [Sjolander96] */
        /* calculated by esl-mixdchlet from 3Di profiles from AlphaFold UniProt for Pfam 35.0, Sean Johnson October 10, 2022 */
  static double defm[9][20] = {
    { 0.1067, 0.0005, 0.1160, 0.0338, 0.1003, 0.0112, 0.0209, 0.0008, 0.0658, 0.0122, 0.0002, 0.0005, 0.0219, 0.0571, 0.0044, 0.0005, 0.0018, 0.0008, 0.0939, 0.0162 },
    { 0.6401, 0.1422, 3.4930, 0.1181, 0.5734, 2.9878, 0.7922, 1.0466, 0.1167, 0.2356, 0.2066, 0.2372, 0.7758, 0.5663, 1.1345, 0.2870, 2.3103, 0.2428, 0.2374, 0.1340 },
    { 0.1697, 0.2209, 0.4218, 0.0144, 0.1101, 0.0217, 0.0470, 0.0011, 0.0040, 0.3621, 0.0002, 0.0970, 0.6368, 0.4338, 0.0770, 0.0979, 0.0027, 0.6304, 0.0249, 0.0029 },
    { 0.0138, 0.0009, 0.1255, 0.0019, 0.0077, 0.1569, 0.0934, 0.1412, 0.0006, 0.0018, 0.0602, 0.0104, 0.0203, 0.0073, 0.0698, 0.0085, 0.1519, 0.0030, 0.0014, 0.0344 },
    { 0.0012, 0.0608, 0.0482, 0.0001, 0.0008, 0.0002, 0.0004, 0.0001, 0.0001, 0.1270, 0.0001, 0.0098, 0.0819, 0.0619, 0.0009, 0.0243, 0.0001, 0.3075, 0.0001, 0.0001 },
    { 0.0870, 0.1654, 0.9139, 0.0081, 0.0526, 0.2830, 0.3382, 0.0202, 0.0029, 0.1051, 0.0014, 0.3509, 0.9256, 0.1516, 0.5037, 0.3548, 0.0918, 0.5718, 0.0182, 0.0403 },
    { 1.4519, 0.0950, 4.6168, 0.0770, 1.2870, 0.2580, 0.2834, 0.0768, 0.3287, 0.1938, 0.0527, 0.1312, 0.6392, 0.9591, 0.2069, 0.1238, 0.1256, 0.1605, 0.7468, 0.0707 },
    { 0.0919, 0.0410, 3.9701, 0.0073, 0.0319, 0.0515, 0.0807, 0.0013, 0.0007, 0.0093, 0.0002, 0.0193, 0.5946, 0.0588, 0.0371, 0.0093, 0.0062, 0.1543, 0.0073, 0.0214 },
    { 0.0013, 0.0646, 0.0522, 0.0004, 0.0007, 0.0105, 0.0147, 0.0006, 0.0003, 0.0417, 0.0005, 0.1561, 0.0977, 0.0235, 0.1144, 0.4455, 0.0008, 0.3284, 0.0005, 0.0006 }
  };

  ESL_ALLOC(pri, sizeof(P7_PRIOR));
  pri->tm = pri->ti = pri->td = pri->em = pri->ei = NULL;

  pri->tm = esl_mixdchlet_Create(1, 3);	 /* single component; 3 params */
  pri->ti = esl_mixdchlet_Create(1, 2);	 /* single component; 2 params */
  pri->td = esl_mixdchlet_Create(1, 2);	 /* single component; 2 params */
  pri->em = esl_mixdchlet_Create(9, 20); /* 9 component; 20 params */
  pri->ei = esl_mixdchlet_Create(1, 20); /* single component; 20 params */

  if (pri->tm == NULL || pri->ti == NULL || pri->td == NULL || pri->em == NULL || pri->ei == NULL) goto ERROR;

  /* Transition priors: originally from Graeme Mitchison. Notes are lost, but we believe
   * they were trained on an early version of Pfam. 
   */
  pri->tm->q[0]        = 1.0;
  pri->tm->alpha[0][0] = 0.7939; /* TMM */
  pri->tm->alpha[0][1] = 0.0278; /* TMI */ /* Markus suggests ~10x MD, ~0.036; test! */
  pri->tm->alpha[0][2] = 0.0135; /* TMD */ /* Markus suggests 0.1x MI, ~0.004; test! */

  pri->ti->q[0]        = 1.0;
  pri->ti->alpha[0][0] = 0.1551; /* TIM */
  pri->ti->alpha[0][1] = 0.1331; /* TII */

  pri->td->q[0]        = 1.0;
  pri->td->alpha[0][0] = 0.9002; /* TDM */
  pri->td->alpha[0][1] = 0.5630; /* TDD */

  /* Match emission priors 
   */  
  for (k = 0; k < 9; k++)
    {
      pri->em->q[k] = defmq[k];
      esl_vec_DCopy(defm[k], 20, pri->em->alpha[k]);
    }

  /* calculated from 3Di profiles from AlphaFold for Pfam 35.0, Sean Johnson June 4, 2022 */
  pri->ei->q[0] = 1.0;
  pri->ei->alpha[0][0]  = 325.;        /* A */
  pri->ei->alpha[0][1]  = 367.;        /* C */
  pri->ei->alpha[0][2]  = 1314.;       /* D */
  pri->ei->alpha[0][3]  = 147.;        /* E */
  pri->ei->alpha[0][4]  = 238.;        /* F */
  pri->ei->alpha[0][5]  = 267.;        /* G */
  pri->ei->alpha[0][6]  = 278.;        /* H */
  pri->ei->alpha[0][7]  = 215.;        /* I */
  pri->ei->alpha[0][8]  = 202.;        /* K */
  pri->ei->alpha[0][9]  = 941.;        /* L */
  pri->ei->alpha[0][10] = 90.;         /* M */
  pri->ei->alpha[0][11] = 310.;        /* N */
  pri->ei->alpha[0][12] = 829.;        /* P */
  pri->ei->alpha[0][13] = 522.;        /* Q */
  pri->ei->alpha[0][14] = 329.;        /* R */
  pri->ei->alpha[0][15] = 890.;        /* S */
  pri->ei->alpha[0][16] = 180.;        /* T */
  pri->ei->alpha[0][17] = 2182.;       /* V */
  pri->ei->alpha[0][18] = 211.;        /* W */
  pri->ei->alpha[0][19] = 163.;        /* Y */

  return pri;

 ERROR:
  if (pri != NULL) p7_prior_Destroy(pri);
  return NULL;
}

/* Function:  p7_prior_CreateNucleic()
 *
 * Purpose:   Creates the default mixture Dirichlet prior for nucleotide
 *            sequences.
 *
 *            The transition priors (match, insert, delete) are all
 *            single Dirichlets, trained on a portion of the rmark dataset
 *
 *            The match emission prior is an eight-component mixture
 *            trained against a portion of the rmark dataset
 *
 *            The insert emission prior is a single Dirichlet with
 *            high $|\alpha|$, such that insert emission probabilities
 *            are essentially fixed by the prior, regardless of
 *            observed count data.
 *
 * Returns:   a pointer to the new <P7_PRIOR> structure.
 */
P7_PRIOR *
p7_prior_CreateNucleic(void)
{
  int status;
  P7_PRIOR *pri = NULL;
  int q;

  /* Plus-1 Laplace prior
  int num_comp = 1;
  static double defmq[2] =  { 1.0  };
  static double defm[1][4] = {
         { 1.0, 1.0, 1.0, 1.0} //
  };
*/

  /* Match emission priors are trained on Rmark3 database
   * Xref: ~wheelert/notebook/2011/0325_nhmmer_new_parameters
   */
   int num_comp = 4;
   static double defmq[4] = { 0.24, 0.26, 0.08, 0.42  };
   static double defm[4][4] = {
         { 0.16,  0.45,  0.12,   0.39},
         { 0.09,  0.03,  0.09,   0.04},
         { 1.29,  0.40,  6.58,   0.51},
         { 1.74,  1.49,  1.57,   1.95}
	};



  ESL_ALLOC(pri, sizeof(P7_PRIOR));
  pri->tm = pri->ti = pri->td = pri->em = pri->ei = NULL;

  pri->tm = esl_mixdchlet_Create(1, 3);        // match transitions; single component; 3 params
  pri->ti = esl_mixdchlet_Create(1, 2);        // insert transitions; single component; 2 params
  pri->td = esl_mixdchlet_Create(1, 2);        // delete transitions; single component; 2 params
  pri->em = esl_mixdchlet_Create(num_comp, 4); // match emissions; X component; 4 params
  pri->ei = esl_mixdchlet_Create(1, 4);        // insert emissions; single component; 4 params

  if (pri->tm == NULL || pri->ti == NULL || pri->td == NULL || pri->em == NULL || pri->ei == NULL) goto ERROR;

  /* Transition priors: roughly, learned from rmark benchmark - hand-beautified (trimming overspecified significant digits)
   */
  pri->tm->q[0]        = 1.0;
  pri->tm->alpha[0][0] = 2.0;  // TMM
  pri->tm->alpha[0][1] = 0.1;  // TMI
  pri->tm->alpha[0][2] = 0.1;  // TMD

  pri->ti->q[0]        = 1.0;
  pri->ti->alpha[0][0] = 0.12; // TIM -  was 0.06 (TW changed 3/19/15)
  pri->ti->alpha[0][1] = 0.4;  // TII -  was 0.2  (TW changed 3/19/15)

  pri->td->q[0]        = 1.0;
  pri->td->alpha[0][0] = 0.5;  // TDM -  was 0.1 (TW changed 3/19/15)
  pri->td->alpha[0][1] = 1.0;  // TDD -  was 0.2 (TW changed 3/19/15)


  /* Match emission priors  */
  for (q = 0; q < num_comp; q++)
    {
      pri->em->q[q] = defmq[q];
      esl_vec_DCopy(defm[q], 4, pri->em->alpha[q]);
    }

  /* Insert emission priors. Should alphas be lower? higher?
   */
  pri->ei->q[0] = 1.0;
  esl_vec_DSet(pri->ei->alpha[0], 4, 1.0);

  return pri;

 ERROR:
  if (pri != NULL) p7_prior_Destroy(pri);
  return NULL;
}



/* Function:  p7_prior_CreateLaplace()
 * Synopsis:  Creates Laplace plus-one prior.
 * Incept:    SRE, Sat Jun 30 09:48:13 2007 [Janelia]
 *
 * Purpose:   Create a Laplace plus-one prior for alphabet <abc>.
 */
P7_PRIOR *
p7_prior_CreateLaplace(const ESL_ALPHABET *abc)
{
  P7_PRIOR *pri = NULL;
  int        status;

  ESL_ALLOC(pri, sizeof(P7_PRIOR));
  pri->tm = pri->ti = pri->td = pri->em = pri->ei = NULL;

  pri->tm = esl_mixdchlet_Create(1, 3);	     /* single component; 3 params */
  pri->ti = esl_mixdchlet_Create(1, 2);	     /* single component; 2 params */
  pri->td = esl_mixdchlet_Create(1, 2);	     /* single component; 2 params */
  pri->em = esl_mixdchlet_Create(1, abc->K); /* single component; K params */
  pri->ei = esl_mixdchlet_Create(1, abc->K); /* single component; K params */

  if (pri->tm == NULL || pri->ti == NULL || pri->td == NULL || pri->em == NULL || pri->ei == NULL) goto ERROR;

  pri->tm->q[0] = 1.0;   esl_vec_DSet(pri->tm->alpha[0], 3,      1.0);  /* match transitions  */
  pri->ti->q[0] = 1.0;   esl_vec_DSet(pri->ti->alpha[0], 2,      1.0);  /* insert transitions */
  pri->td->q[0] = 1.0;   esl_vec_DSet(pri->td->alpha[0], 2,      1.0);  /* delete transitions */
  pri->em->q[0] = 1.0;   esl_vec_DSet(pri->em->alpha[0], abc->K, 1.0);  /* match emissions    */
  pri->ei->q[0] = 1.0;   esl_vec_DSet(pri->ei->alpha[0], abc->K, 1.0);  /* insert emissions   */
  return pri;

 ERROR:
  p7_prior_Destroy(pri);
  return NULL;
}


/* Function:  p7_prior_Destroy()
 * Incept:    SRE, Sat Mar 24 09:55:09 2007 [Janelia]
 *
 * Purpose:   Frees a mixture Dirichlet prior.
 */
void
p7_prior_Destroy(P7_PRIOR *pri)
{
  if (pri == NULL) return;
  if (pri->tm != NULL) esl_mixdchlet_Destroy(pri->tm);
  if (pri->ti != NULL) esl_mixdchlet_Destroy(pri->ti);
  if (pri->td != NULL) esl_mixdchlet_Destroy(pri->td);
  if (pri->em != NULL) esl_mixdchlet_Destroy(pri->em);
  if (pri->ei != NULL) esl_mixdchlet_Destroy(pri->ei);
  free(pri);
}



/* Function:  p7_ParameterEstimation()
 * Incept:    SRE, Sat Mar 24 10:15:37 2007 [Janelia]
 *
 * Purpose:   Given an <hmm> containing weighted counts, and
 *            a mixture Dirichlet prior <pri>: calculate mean
 *            posterior parameter estimates for all model parameters,
 *            converting the HMM to a parameterized probabilistic
 *            model.
 *            
 *            If <pri> is <NULL>, then model parameters are calculated
 *            as frequencies, by normalization of <hmm>.
 *            
 * Args:      hmm - profile structure, containing counts.
 *            pri - mixture Dirichlet prior structure, or <NULL>.
 *            
 * Returns:   <eslOK> on success.
 */
int
p7_ParameterEstimation(P7_HMM *hmm, const P7_PRIOR *pri)
{
  int   k;
  double c[p7_MAXABET];
  double p[p7_MAXABET];

  //  double mix[p7_MAXDCHLET];
  
  /* Special case of pri=NULL: convert to frequencies.*/  
  if (pri==NULL) return p7_hmm_Renormalize(hmm);

  /* Match transitions 0,1..M: 0 is the B state
   * TMD at node M is 0.
   */
  for (k = 0; k <= hmm->M; k++) {
    esl_vec_F2D(hmm->t[k], 3, c);
    esl_mixdchlet_MPParameters(pri->tm, c, p);
    esl_vec_D2F(p, 3, hmm->t[k]);
  }
  hmm->t[hmm->M][p7H_MD] = 0.0;
  esl_vec_FNorm(hmm->t[hmm->M], 3);

  /* Insert transitions, 0..M
   */
  for (k = 0; k <= hmm->M; k++) {
    esl_vec_F2D(hmm->t[k]+3, 2, c);
    esl_mixdchlet_MPParameters(pri->ti, c, p);
    esl_vec_D2F(p, 2, hmm->t[k]+3);
  }

  /* Delete transitions, 1..M-1
   * For k=0, which is unused; convention sets TMM=1.0, TMD=0.0
   * For k=M, TMM = 1.0 (to the E state) and TMD=0.0 (no next D; must go to E).
   */
  for (k = 1; k < hmm->M; k++) {
    esl_vec_F2D(hmm->t[k]+5, 2, c);
    esl_mixdchlet_MPParameters(pri->td, c, p);
    esl_vec_D2F(p, 2, hmm->t[k]+5);
  }
  hmm->t[0][p7H_DM] = hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[0][p7H_DD] = hmm->t[hmm->M][p7H_DD] = 0.0;

  /* Match emissions, 1..M
   * Convention sets mat[0] to a valid pvector: first elem 1, the rest 0.
   */
  for (k = 1; k <= hmm->M; k++) {
    esl_vec_F2D(hmm->mat[k], hmm->abc->K, c);
    esl_mixdchlet_MPParameters(pri->em, c, p);
    esl_vec_D2F(p, hmm->abc->K, hmm->mat[k]);
  }
  esl_vec_FSet(hmm->mat[0], hmm->abc->K, 0.);
  hmm->mat[0][0] = 1.0;

  /* Insert emissions 0..M
   */
  for (k = 0; k <= hmm->M; k++) {
    esl_vec_F2D(hmm->ins[k], hmm->abc->K, c);
    esl_mixdchlet_MPParameters(pri->ei, c, p);
    esl_vec_D2F(p, hmm->abc->K, hmm->ins[k]);
  }
  return eslOK;
}


