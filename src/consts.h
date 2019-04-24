#ifndef _CPPTEST_CONSTS_H_
#define _CPPTEST_CONSTS_H_

#define SIM_YEARS

#define AGE_START 15
#define DISEASE_STATUS 2 // pDS - HIV positive, HIV negative
#define SEXES 2  // NG - Male, Female
#define MODEL_AGES 66 // pAG - 15, 16, ..., 78, 79, 80, 80+
#define PROJECTION_YEARS 56 // PROJ_YEARS - 1970, 1971, ...
#define AGE_GROUPS 9 // hAG 15-19, 20-24, 25-29, ..., 45-49, 50+
#define CD4_STAGES 7 // hDS CD4 count stages
#define TREATMENT_STAGES 3// hTS No of treatment stages including untreated

#define FERT_AGES 35 // pAG_FERT
#define pIDX_FERT 0 // pIDX_FERT

#define FERT_AGE_GROUPS 8 // hAG_FERT No of age groups for fertility
#define hIDX_FERT 0 // hIDX_FERT 

// Below constants used for accessing data
#define MALE 0
#define FEMALE 1

#define HIVN 0
#define HIVP 1

#endif