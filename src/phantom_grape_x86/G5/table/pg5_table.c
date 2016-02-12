#include<math.h>
#include<stdio.h>
#include "sse_type.h"
#define GLOBAL_DEFINE
#include "pg5_table.h"
#undef GLOBAL_DEFINE

static double reff(double sft_length, double radius)
{
	double ret;
	double eta  = 2.0*radius/sft_length;

	if(eta < 1.0){
		ret = (224.+eta*eta*(-224.+eta*(70.+eta*(48.-eta*21.)))) 
			* 2.e0/(35.0*pow(sft_length,3));
	}else if(eta < 2.0){
		ret = (12./(eta*eta)-224.+eta*(896.+eta*(-840.+eta*(224.+eta*(70.+eta*(-48.+eta*7.)))))) *2.e0/(35.0*eta*pow(sft_length,3));
	}else{
		ret = 1.0/pow(radius,3);
	}
	return ret;
}

static double Sft_for_PP = SFT_FOR_PP;
static double Sft_for_PM = SFT_FOR_PM;

static double s2_force(double r){
		return (reff(Sft_for_PP, r)  - reff(Sft_for_PM, r));
}

void pg5_gen_s2_force_table(double sft_for_PP, double sft_for_PM){
	Sft_for_PP = sft_for_PP;
	Sft_for_PM = sft_for_PM;
	pg5_gen_force_table(s2_force, sft_for_PM);
}

static float R2scacle;
union pack32{
	float f;
	unsigned u;
};

void pg5_gen_force_table(
	double (*force_func)(double r), // caluclate r^{-3} from r in N-body unit
	double rcut // cut-off radius in N-body unit
) {
	int i;
	const unsigned tick = (1<<(23-FRC_BIT));
	union pack32 m32;
	const float r2max = (rcut * rcut);
	// const float r2scale = ((1<<(1+(1<<EXP_BIT))) - 3) / r2max;
	// const float r2scale = ((1<<(1+(1<<EXP_BIT))) - 3 - (1<<FRC_BIT)) / r2max;
	const float fmax = (1<<(1<<EXP_BIT)) * (2.0 - 1.0/(1<<FRC_BIT));
	const float r2scale = (fmax - 2.0f) / r2max;
		// so that 0 <= r2/r2max <= 2^17 - 3
		// or 2 <= s < 2^17, where s = r2/r2max + 2.0
		// in binary, 0_100000000_00000000... <= s <= 0_100001111_11110000...

	R2scacle = r2scale;
#ifndef DEBUG_MODE
	pg5_set_xscale((double)sqrtf(r2scale));
#endif
	// table value and...
	for(i=0, m32.f = 2.0f; 
		i<TBL_SIZE; 
		i++, m32.u += tick)
	{
		float f = m32.f;
		float r2 = (f - 2.0) / r2scale;
		float r = sqrtf(r2);
		Force_table[i][0] = force_func(r);
		//#if 0000
#if 1
		printf("%+e %+e\n", r/SFT_FOR_PM, Force_table[i][0] * SFT_FOR_PM * SFT_FOR_PM * SFT_FOR_PM);
#endif
	}
#if 0
	// BIAS correcton
	for(i=0, m32.f = 2.0f; 
		i<TBL_SIZE; 
		i++)
	{
		static float delta_backward;
		float x0 = m32.f;
		m32.u += tick;
		float x1 = m32.f;
		float f = (x0+x1)/2.0;
		float r2 = (f - 2.0) / r2scale;
		float r = sqrtf(r2);
		float delta_forward = 
			force_func(r) - (Force_table[i][0] + Force_table[i+1][0])/2.0;
		if (i==0) delta_backward = delta_forward;
		float delta = (delta_forward + delta_backward)/2.0;
		if (i==TBL_SIZE-1) delta = 0.0;

		Force_table[i][0] += 2./3. * delta;
		delta_backward = delta_forward;
	}
#endif
	// ...slope
	for(i=0, m32.f = 2.0f; 
		i<TBL_SIZE-1; 
		i++)
	{
		float x0 = m32.f;
		m32.u += tick;
		float x1 = m32.f;
		float y0 = Force_table[i][0];
		float y1 = (i==TBL_SIZE-1) ? 0.0 : Force_table[i+1][0];
		Force_table[i][1] = (y1-y0)/(x1-x0);
	}
	Force_table[i][1] = 0.0f;
}

// The followings are for debug
#ifdef DEBUG_MODE
#include <stdio.h>
#include <assert.h>
static void dump_table(FILE *fp){
	int i;
	const unsigned tick = (1<<(23-FRC_BIT));
	union pack32 m32;
	float r2scale = R2scacle;

	for(i=0, m32.f = 2.0f; 
		i<TBL_SIZE; 
		i++, m32.u += tick)
	{
		float f = m32.f;
		float r2 = (f - 2.0) / r2scale;
		float r = sqrtf(r2);
		float force = r * Force_table[i][0];
		fprintf(fp, "%e %e %e %e\n", r, force, Force_table[i][0], Force_table[i][1]);
	}
}

static float refer_table(float r){
	union pack32 m32;
	float x = r*r*R2scacle + 2.0f;
	m32.f = x;
	int idx = (m32.u >> (23-FRC_BIT)) & (TBL_SIZE-1);
	assert(idx < TBL_SIZE);
	m32.u >>= (23-FRC_BIT);
	m32.u <<= (23-FRC_BIT);
	float dx = (x - m32.f);
	return Force_table[idx][0] + dx*Force_table[idx][1];
	// return Force_table[idx][0];
}	

#define NUM_POINT 10000
int main(){
	int i;
	float tick = SFT_FOR_PM/NUM_POINT;

	FILE *fp1 = fopen("force_table.dat" ,"w");
	FILE *fp2 = fopen("force_out.dat" ,"w");

	pg5_gen_force_table(s2_force, SFT_FOR_PM);
	dump_table(fp1);

	for(i=0; i<=NUM_POINT; i++){
		float r = i*tick;
		double force = r * refer_table(r);
		double realf = r * s2_force(r);
		double ferror = (force - realf);
		fprintf(fp2, "%e %e %e %e %e\n", r, force, realf, ferror, ferror/force);
	}
	fclose(fp1);
	fclose(fp2);
	return 0;
}
#endif
