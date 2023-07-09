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
        /* calculated by esl-mixdchlet from 3Di profiles from AlphaFold UniProt for Pfam 35.0, Sean Johnson October 10, 2022 */
  static double defmq[9] = {  0.0674, 0.0670, 0.0696, 0.1586, 0.0852, 0.0746, 0.1386, 0.0406, 0.2984 };

				/* default match mixture Dirichlet components [Sjolander96] */
        /* calculated by esl-mixdchlet from 3Di profiles from AlphaFold UniProt for Pfam 35.0, Sean Johnson October 10, 2022 */
  static double defm[9][20] = {
    { 0.0717, 0.0008, 0.0812, 0.0327, 0.0959, 0.0041, 0.0034, 0.0006, 0.1050, 0.0023, 0.0005, 0.0005, 0.0075, 0.0252, 0.0009, 0.0006, 0.0007, 0.0012, 0.1320, 0.0012 },
    { 0.0374, 0.0007, 0.1979, 0.0144, 0.0120, 0.1212, 0.1412, 0.0027, 0.0004, 0.0018, 0.0001, 0.0075, 0.0185, 0.0067, 0.0519, 0.0046, 0.0403, 0.0014, 0.0013, 0.0751 },
    { 0.0470, 0.0665, 0.4273, 0.0042, 0.0288, 0.3072, 0.2869, 0.0410, 0.0039, 0.0473, 0.0060, 0.2945, 0.3659, 0.0720, 0.5856, 0.3434, 0.1153, 0.2149, 0.0100, 0.0331 },
    { 0.0222, 0.2322, 0.2737, 0.0016, 0.0132, 0.0104, 0.0235, 0.0014, 0.0005, 0.1942, 0.0002, 0.1189, 0.4826, 0.1504, 0.0795, 0.1709, 0.0027, 1.5896, 0.0020, 0.0021 },
    { 0.2604, 0.3586, 3.0590, 0.0237, 0.1460, 0.1908, 0.2672, 0.0210, 0.0143, 0.3033, 0.0083, 0.2853, 3.2179, 0.4208, 0.3264, 0.2880, 0.0624, 1.3300, 0.0573, 0.0290 },
    { 0.4267, 0.0613, 0.4502, 0.0412, 0.2849, 0.0462, 0.0600, 0.0038, 0.0225, 0.1960, 0.0014, 0.0317, 0.2617, 0.5496, 0.0480, 0.0259, 0.0120, 0.1074, 0.0892, 0.0091 },
    { 0.0633, 0.0353, 4.5152, 0.0047, 0.0208, 0.0317, 0.0496, 0.0025, 0.0016, 0.0104, 0.0012, 0.0157, 0.7239, 0.0386, 0.0259, 0.0105, 0.0077, 0.1595, 0.0058, 0.0156 },
    { 0.0053, 0.0007, 0.0977, 0.0008, 0.0039, 0.1103, 0.0249, 0.5747, 0.0013, 0.0007, 0.1872, 0.0013, 0.0113, 0.0038, 0.0305, 0.0022, 0.3149, 0.0013, 0.0019, 0.0023 },
    { 0.0012, 0.0374, 0.0221, 0.0001, 0.0010, 0.0005, 0.0010, 0.0001, 0.0001, 0.0812, 0.0001, 0.0295, 0.0467, 0.0389, 0.0169, 0.0749, 0.0001, 0.2099, 0.0001, 0.0001 }
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

  /* calculated from 3Di profiles from AlphaFold for Pfam 35.0, Sean Johnson October 10, 2022 */
  pri->ei->q[0] = 1.0;
  pri->ei->alpha[0][0]  = 322.;         /* A */
  pri->ei->alpha[0][1]  = 365.;         /* C */
  pri->ei->alpha[0][2]  = 1357.;        /* D */
  pri->ei->alpha[0][3]  = 145.;         /* E */
  pri->ei->alpha[0][4]  = 236.;         /* F */
  pri->ei->alpha[0][5]  = 265.;         /* G */
  pri->ei->alpha[0][6]  = 276.;         /* H */
  pri->ei->alpha[0][7]  = 213.;         /* I */
  pri->ei->alpha[0][8]  = 199.;         /* K */
  pri->ei->alpha[0][9]  = 930.;         /* L */
  pri->ei->alpha[0][10] = 89.;          /* M */
  pri->ei->alpha[0][11] = 307.;         /* N */
  pri->ei->alpha[0][12] = 845.;         /* P */
  pri->ei->alpha[0][13] = 517.;         /* Q */
  pri->ei->alpha[0][14] = 327.;         /* R */
  pri->ei->alpha[0][15] = 879.;         /* S */
  pri->ei->alpha[0][16] = 178.;         /* T */
  pri->ei->alpha[0][17] = 2182.;        /* V */
  pri->ei->alpha[0][18] = 210.;         /* W */
  pri->ei->alpha[0][19] = 160.;         /* Y */

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


