#include "../measure.hpp"

#include <cstring>

#include "emmintrin.h"
#include "pmmintrin.h"

class FullMeasure_CP_Omega1 : public FullMeasure {
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
    double den = 0;
    int count = 0;

    for (size_t px = 0; px < numX; px++){
    	pc = prior.cdata();
    	pxc = cond.cdata() + count;
    	for (size_t pc_iter = 0; pc_iter < numC ; pc_iter++){
    		num 	+=	*pc * *pc *  *pxc * *pxc;
    		den 	+=	*pc * *pxc;
    		pc++;
    		pxc += numX;
    		if (numX%2)	pxc ++;
    	}
    	res += num/den;
    	count++;
    	num = 0;	den = 0;
    }

    double* result = out.data(); 
    result  = result + pos;
    *result = res;
    return 1;
  }

};

DEFINE_DEFAULT_FACTORY(FullMeasure, FullMeasure_CP_Omega1)
EXPORT_FACTORY(Factory)


