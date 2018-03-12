#include "../measure.hpp"
#include <math.h>
#include <cstring>

#include "emmintrin.h"
#include "pmmintrin.h"

class FullMeasure_CP_Omega3 : public FullMeasure {
public:
  virtual unsigned measure(FullDistribution const& dist, std::vector<double>& out, unsigned pos) const override {
    Prior const& prior = dist.first;
    Conditional const& cond = dist.second;
    Index numC = cond.dimSize(1);
    Index numX = cond.dimSize(0);

    double const* pc = prior.cdata();
    double const* pc_end = pc + numC;
    double const* pxc = cond.cdata();
    double term1 = 0. , term2 = 0. ;
    double num = 0. ,den = 0. , division = 0. ;
	double res = 0. ;
	
// Calculating Term1  =  sqrt(sum_{pr^2(c)})
    while(pc < pc_end){
		term1 += (*pc)*(*pc);
		pc++;
	}
	term1 = sqrt(term1);

// Calculating Term2 
	// used pc_iter so as not to confuse with pc pointer
	for (size_t pc_iter = 0; pc_iter < numC ; pc_iter++){
	    for (size_t px = 0; px < numX; px++){
	    	//pxc 	= pxc + px;
	    	num 	= *pxc * *pxc;
			den = 0;
	    	double const* pxc1 = cond.cdata();
	    	pxc1 	= pxc1 + px ;
	    	for (size_t px1 = 0; px1 < numC; px1++){
	    		den 	+=  *pxc1 * *pxc1 ;
	    		pxc1 	+= numX;
	    		if (numX%2) pxc1++;    		
	    	}
	    	den 	= sqrt(den);
	    	division += num/den;
			pxc++;
	    }
		if (pc_iter == 0	) 	term2 = division;		
		if (term2 > division) 	term2 = division;
	//	pxc ++;
		if(numX % 2) pxc++;
		division = 0.;
	}

	res = term1 * term2;
    double* result = out.data();
    result  = result + pos;
    *result = res;
    return 1;

  }

};

DEFINE_DEFAULT_FACTORY(FullMeasure, FullMeasure_CP_Omega3)
EXPORT_FACTORY(Factory)


