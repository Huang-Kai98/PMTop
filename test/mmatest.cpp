// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-16.
//
#include <GCMMASolver.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

double Squared(double x) { return x*x; }

struct Problem {
	int n, m;
	std::vector<double> x0, xmin, xmax;

	Problem()
		: n(2)
		, m(1)
		, x0({5, 5, 5})
		, xmin(n, -10.0)
		, xmax(n, 10.0)
	{ }

	void Obj(const double *x, double *f0x, double *fx) {
		f0x[0] = x[0]*x[0]+x[1]*x[1];
        fx[0]=x[0]+x[1]-4;
	}

	void ObjSens(const double *x, double *f0x, double *fx, double *df0dx, double *dfdx) {
		Obj(x, f0x, fx);
		df0dx[0]=2*x[0];df0dx[1]=2*x[1];
        dfdx[0]=1;dfdx[1]=1;
	}
};

void Print(double *x, int n, const std::string &name = "x") {
	std::cout << name << ":";
	for (int i=0;i<n;i++) {
		std::cout << " " << x[i];
	}
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {
	std::cout << "///////////////////////////////////////////////////" << std::endl;
	std::cout << "// Test the GCMMA Algorithm" << std::endl;
	std::cout << "///////////////////////////////////////////////////" << std::endl;

	Problem toy;


	double f, fnew;
	std::vector<double> df(toy.n);
	std::vector<double> g(toy.m), gnew(toy.m);
	std::vector<double> dg(toy.n * toy.m);

	std::vector<double> x = toy.x0;
	std::vector<double> xold = x;
	std::vector<double> xnew(toy.n);

	// Print initial values
	toy.Obj(x.data(), &f, g.data());
	std::cout << "f: " << f << std::endl;
	Print(g.data(), toy.m, "g");

	// Initialize GCMMA
	GCMMASolver gcmma(toy.n, toy.m, 0, 1000, 1);
	double ch = 1.0;
	int maxoutit = 8;
	for (int iter = 0; ch > 0.0002 && iter < maxoutit; ++iter) {
		toy.ObjSens(x.data(), &f, g.data(), df.data(), dg.data());
			// GCMMA version
		gcmma.OuterUpdate(xnew.data(), x.data(), f, df.data(),
				g.data(), dg.data(), toy.xmin.data(), toy.xmax.data());

			// Check conservativity
		toy.Obj(xnew.data(), &fnew, gnew.data());
		bool conserv = gcmma.ConCheck(fnew, gnew.data());
			//std::cout << conserv << std::endl;
		for (int inneriter = 0; !conserv && inneriter < 15; ++inneriter) {
				// Inner iteration update
		gcmma.InnerUpdate(xnew.data(), fnew, gnew.data(), x.data(), f,
					df.data(), g.data(), dg.data(), toy.xmin.data(), toy.xmax.data());

				// Check conservativity
		toy.Obj(xnew.data(), &fnew, gnew.data());
		conserv = gcmma.ConCheck(fnew, gnew.data());
				//std::cout << conserv << std::endl;
		}
		x = xnew;


		// Compute infnorm on design change
		ch = 0.0;
		for (int i=0; i < toy.n; ++i) {
			ch = std::max(ch, std::abs(x[i] - xold[i]));
			xold[i] = x[i];
		}

		// Print to screen
		printf("it.: %d, obj.: %f, ch.: %f \n", iter, f, ch);
		Print(x.data(), toy.n);
		toy.Obj(x.data(), &f, g.data());
		std::cout << "f: " << f << std::endl;
		Print(g.data(), toy.m, "g");
		std::cout << std::endl;
	}

	return 0;
}
