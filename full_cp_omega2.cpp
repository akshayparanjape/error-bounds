#include "../measure.hpp"
#include <math.h>
#include <cstring>

#include "emmintrin.h"
#include "pmmintrin.h"

class FullMeasure_CP_Omega2 : public FullMeasure {
public:
  virtual unsigned measure(FullDistribution const& dist, std::vector<double>& out, unsigned pos) const override {
    Prior const& prior = dist.first;
    Conditional const& cond = dist.second;
    Index numC = cond.dimSize(1);
    Index numX = cond.dimSize(0);
    std::vector<double> px_vec(numX);

    double const* pc = prior.cdata();
    double const* pc_end = pc + numC;
    double const* pxc = cond.cdata();

    double res = 0;
    double num = 0;
    double den1 = 0;
    double den2 = 0;
    int count = 0;
    double division = 0;

    for (size_t pc_iter = 0; pc_iter < numC ; pc_iter++){
		den1 += (*pc)*(*pc);
		pc++;
	}
	den1 = sqrt(den1);

	pc = prior.cdata();
	for (size_t pc_iter = 0; pc_iter < numC ; pc_iter++){
	    for (size_t px = 0; px < numX; px++){
	    	num 	= *pxc * *pxc;
			den2 = 0;
	    	double const* pxc1 = cond.cdata();
	    	pxc1 	= pxc1 + px ;	    	
	    	for (size_t px1 = 0; px1 < numC; px1++){
	    		den2 	+=  *pxc1 * *pxc1 ;
	    		pxc1 	+= numX;
	    		if (numX%2) pxc1++;    		
	    	}
	    	den2 	= sqrt(den2);
	    	division += num/den2;
	    	num = 0;
			pxc++;
	    }
	    res += (*pc) * (*pc) * division;
		division = 0.;
		pc++;
		if(numX % 2) pxc++;
	}	

	res = res/den1;
    double* result = out.data(); 
    result  = result + pos;
    *result = res;
    return 1;
  }

};

DEFINE_DEFAULT_FACTORY(FullMeasure, FullMeasure_CP_Omega2)
EXPORT_FACTORY(Factory)


