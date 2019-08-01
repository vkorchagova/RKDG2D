#include "Flux.h"

#include <algorithm>

using namespace std;

double Flux::c_av(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const
{
	double semiRho = 0.5*(solOne[0] + solTwo[0]);
	double semiP = 0.5*(phs.getPressure(solOne) + phs.getPressure(solTwo));

	return sqrt(phs.cpcv * semiP / semiRho);
} // end c for edge

numvector<double, dimPh> Flux::lambdaF_Roe(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const
{
	double sqrtRhoLeft = sqrt(solOne[0]);
	double sqrtRhoRight = sqrt(solTwo[0]);

	double sumSqrtRho = sqrtRhoLeft + sqrtRhoRight;

	double u_av = (solOne[1] / sqrtRhoLeft + solTwo[1] / sqrtRhoRight) / sumSqrtRho;
	double h_av = ((solOne[4] + phs.getPressure(solOne)) / sqrtRhoLeft + (solTwo[4] + phs.getPressure(solTwo)) / sqrtRhoRight) / sumSqrtRho;

	double c_av = sqrt((phs.cpcv - 1) * (h_av - 0.5 * sqr(u_av)));

	//    if (getPressure(solOne) < 0 || getPressure(solTwo) < 0)
	//        cout << "kkk";

	return{ u_av - c_av, u_av, u_av, u_av, u_av + c_av };

}

numvector<double, dimPh> Flux::lambdaF_Einfeldt(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const
{
	double cLeft = phs.c(solOne);
	double cRight = phs.c(solTwo);

	double uLeft = solOne[1] / solOne[0];
	double uRight = solTwo[1] / solTwo[0];

	double sqrtRhoLeft = sqrt(solOne[0]);
	double sqrtRhoRight = sqrt(solTwo[0]);
	double sumSqrtRho = sqrtRhoLeft + sqrtRhoRight;

	double eta2 = 0.5 * sqrtRhoLeft * sqrtRhoRight / sqr(sqrtRhoLeft + sqrtRhoRight);
	double u_av = (solOne[1] / sqrtRhoLeft + solTwo[1] / sqrtRhoRight) / sumSqrtRho;
	double c_av = sqrt((sqrtRhoLeft * sqr(cLeft) + sqrtRhoRight * sqr(cRight)) / sumSqrtRho + \
		eta2 * sqr(uRight - uLeft));

	return{ u_av - c_av, u_av, u_av, u_av, u_av + c_av };

}

numvector<double, dimPh> Flux::lambdaF_Toro(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const
{
	double uLeft = solOne[1] / solOne[0];
	double uRight = solTwo[1] / solTwo[0];

	double cLeft = phs.c(solOne);
	double cRight = phs.c(solTwo);

	double pLeft = phs.getPressure(solOne);
	double pRight = phs.getPressure(solTwo);

	double pvrs = 0.5 * (pLeft + pRight) + \
		0.125 * (uLeft - uRight) * (solOne[0] + solTwo[0]) * (cLeft + cRight);

	double pStar = max(0.0, pvrs);

	double qLeft = pStar > pLeft ? \
		sqrt(1.0 + 0.5 * (phs.cpcv + 1.0) * (pStar / pLeft - 1.0) / phs.cpcv) : \
		1.0;

	double qRight = pStar > pRight ? \
		sqrt(1.0 + 0.5 * (phs.cpcv + 1.0) * (pStar / pRight - 1.0) / phs.cpcv) : \
		1.0;

	return{ uLeft - qLeft * cLeft, uLeft, uLeft, uLeft, uRight + qRight * cRight };
}

numvector<double, dimPh> Flux::lambdaF_semisum(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const
{
	double u = 0.5*(solOne[1] / solOne[0] + solTwo[1] / solTwo[0]);
	double soundSpeed = c_av(solOne, solTwo);

	return{ u - soundSpeed, u, u, u, u + soundSpeed };
} // end lambdaF

numvector<double, dimPh> Flux::lambdaF(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const
{
	return lambdaF_Roe(solOne, solTwo);
} // end lambdaF

